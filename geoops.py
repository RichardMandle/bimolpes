# functions for geometry operations (geoops)

import numpy as np
import os 
from scipy.spatial.transform import Rotation as R
from rdkit import Chem
import glob

# grab our own modules
import gauops as gps
import orcops as ops
import processing as pro

def get_geometries(path):
    """
    NEW - a general geometry loader.

    based on the extension, calls the right file handling function

    Returns:
        List of geometries (each geometry is a list of coordinate lines)
    """

    if path is None:
        raise ValueError("No geometry file specified.")

    # need to handle glob patterns explicitly
    if any(ch in path for ch in ["*", "?", "["]):
        # assume XYZ glob for now
        # probably fair enough
        return ops.get_xyz_geometries(path)

    ext = os.path.splitext(path)[1].lower()

    if ext in [".log", ".out"]:
        return gps.get_gau_geometries(path)

    elif ext == ".xyz":
        return ops.get_xyz_geometries(path)

    # aim here is we could support all sorts of files quite easily.

    else:
        raise ValueError(f"Unsupported geometry format: {ext}")

def split_to_xyz(values_to_split):
    '''
    Little function that splits a list/array of 1-2-3 values into something we can parse later on (either res. or disp. data?)
    '''
    split_values = [np.abs(values_to_split[0]), np.abs(values_to_split[0]), np.abs(values_to_split[0])]  # if 1 passed, x=y=z; take the absolute incase the user passes a negative value
    if len(values_to_split) > 1:
        split_values[1] = split_values[2] = values_to_split[1]                   # if 2 passed, x!=y=z
    if len(values_to_split) > 1:
        split_values[2] = values_to_split[2]                                     # if 3 res passed, x!=y!=z
        
    return np.array(split_values)

def write_minima_files(args, minima_list):
    """
    Write Gaussian input files for each minima in the minima_list.
    
    Args:
        args: Parsed command line arguments.
        minima_list: List of minima with their coordinates and energies.
    """
    if args.mol is None:
        raise ValueError("Error: Ye have to specify a primary molecule file (-mol)")
        
    #get an output path to write to
    file_name = f'{args.filename}_min' # append _min sow we know we are dealing with some generated minima
    os.makedirs(file_name, exist_ok=True)
    pro.clean_path(file_name) # clean that directory
    
    for i, minima in enumerate(minima_list):
        count = i + 1
        displacement = f'dx={minima[0]}/dy={minima[1]}/dz={minima[2]}/rx={minima[3]}/ry={minima[4]}/rz={minima[5]}/'
        this_file_name = os.path.join(file_name, os.path.basename(f'{file_name}_{count}'))

        # get geometries from user supplied file(s)
        frag1 = get_geometries(args.mol)[-1]
        frag2 = get_geometries(args.mol2)[-1] if args.mol2 else frag1
        
        # Look at rotation in the displacement line of our .log file; replicate it here
        if (np.abs(minima[3]) + np.abs(minima[4]) + np.abs(minima[5])) != 0:
            frag2 = rotate_coordinates(args, frag2,  x_angle=minima[3], y_angle=minima[4], z_angle=minima[5])

        if args.orca:
            # if we're using ORCA, call the ORCA specific routines
            ops.write_inp(args, frag1, frag2, this_file_name)                                                                # TO DO - function not written.
            ops.make_sge_job(args, outname = os.path.join(file_name, os.path.basename(file_name)), startjob=1, endjob=count) # make the SGE job for ORCA (needs testing)
            print(f'Wrote {count} minima to new ORCA .inp files in {file_name}')    
            
        # add fragment labels 
        if not args.orca:
            # if we aren't using ORCA, do it this way:
            frag1, atom_count = add_fragment_label(args, frag1, 1, atom_count = 0)
            frag2, atom_count = add_fragment_label(args, translate_coordinates(frag2, dx=minima[0], dy=minima[1], dz=minima[2]), 2, atom_count = atom_count)
            
            gps.write_gjf(args, frag1, frag2, displacement, this_file_name)
            gps.make_sge_job(args, outname = os.path.join(file_name, os.path.basename(file_name)), startjob=1, endjob=count) # make the SGE job
            print(f'Wrote {count} minima to new Gaussian .gjf files in {file_name}')
            

    
def write_grid(args):
    '''
    Master function for grid writing.
    
    Generates a grid of displacements, handling calls to create_inputs
    and writing to .gjf (or .inp), .sh and .zip files as needed.
    '''
    outpath = pro.get_output_filename(args.out, args.mol, args.mol2)
    outpath = os.path.abspath(outpath)
       
    if os.path.exists(outpath):
        if not args.backoff:
            pro.backup_existing_directory(outpath)
        if args.backoff:
            print(f'WARNING: backup of {outpath} has been overridden by -backoff flag. Overwriting...')
    
    os.makedirs(outpath, exist_ok=True)
    file_plus_path = os.path.join(outpath, os.path.basename(outpath))

    if args.est:
        args.x, args.y, args.z = estimate_grid_from_mol(args)          # here, estimate the x/y/z ranges from the input file (args.mol) using estimate_grid_from_mol.
        
    grid = gen_grid(args)     # get a grid using the user supplied spacings.
    
    outname, count =  create_inputs(args, grid = grid, outname = file_plus_path)
    if args.frz:
        print(f'Freezing Atoms: {args.frz}')
        
    if args.orca:    
        ops.make_sge_job(args, outname = file_plus_path, startjob=1, endjob=count) # make the SGE job
           
    if not args.orca:    
        gps.make_sge_job(args, outname = file_plus_path, startjob=1, endjob=count) # make the SGE job

    if args.zip: # bundle into .zip for unleashing on HPC (quicker to upload one file than 100k small ones)
        pro.zip_files(args, dir_path = outpath, zip_file = outpath)

# TO DO - would be good to simplify this function and just read the arguments directly.
import numpy as np

def create_inputs(args, grid, count=1, outname=None):
    ''' 
    Generates translated molecular geometries and filters based on distance constraints.
    Optimized using NumPy vectorization, but keeps molecule data in text form.
    '''
    
    if outname is None:
        outname = os.path.splitext(os.path.basename(args.mol))[0] + '_grid'
    
    pro.clean_path(outname)  # Clean directory before writing
    
    too_close = too_far = 0  # Counters for rejected structures

    # Load geometries only once
    geometry = get_geometries(args.mol)[-1]
    frag2_geometry = geometry  # Assume identical unless mol2 is provided

    if args.mol2 and args.mol2 != args.mol:
        frag2_geometry = get_geometries(args.mol2)[-1]

    # Apply rotation once, outside the loop
    rot_xyz = split_to_xyz(args.rot)
    if np.any(rot_xyz):  # If any rotation is nonzero
        frag2_geometry = rotate_coordinates(args, frag2_geometry, *rot_xyz)

    frag1, atom_count = add_fragment_label(args, geometry, 1, atom_count=0)
    
    # **Precompute translations** (avoid calling `translate_coordinates` inside the loop)
    translated_geometries = {
        tuple(gri): translate_coordinates(frag2_geometry, dx=gri[0], dy=gri[1], dz=gri[2])
        for gri in grid
    }

    for gri, frag2 in translated_geometries.items():
        # **Vectorized distance checks (works directly with text-based fragments)**
        if not check_fragments_too_close(frag1, frag2, min_cutoff=args.min_dist):
            if not check_fragments_too_far(frag1, frag2, max_cutoff=args.max_dist):
                
                displacement = f'dx={gri[0]}/dy={gri[1]}/dz={gri[2]}/rx={rot_xyz[0]}/ry={rot_xyz[1]}/rz={rot_xyz[2]}'
                file_name = f"{outname}_{count}"

                if args.orca:
                    ops.write_inp(args, geometry, frag2, displacement=displacement, file_name=file_name)
                else:
                    frag2, atom_count = add_fragment_label(args, frag2, 2, atom_count = atom_count) # add fragment #2 labels
                    gps.write_gjf(args, frag1, frag2, displacement=displacement, file_name=file_name)
                
                count += 1
            else:
                too_far += 1
        else:
            too_close += 1

    print(f'{too_close} geometries rejected for atom-atom close contacts <= {args.min_dist} Angstrom')
    print(f'{too_far} geometries rejected for minimum atom-atom distance > {args.max_dist} Angstrom')
    print(f"A total of {count-1} geometries were written to {'ORCA .inp format' if args.orca else 'Gaussian .gjf format'}")

    return outname, count


def estimate_grid_from_mol(args):
    '''
    Reads incoming geometry (mol via args.mol) and estimate 
    a reasonable set of  x/y/z limits as the floor/ceiling of the min/max x/y/z coordinates of the input geometry.
    
    Args:
        args        -   arguments passed by bimolpes
        both_hemispheres (bool): If False, truncate at z = 0. If True, do both +/-z and get a 'shell' of complexation energies.
        dx/dy/dz: additional user requested displacement in X/Y/Z, respectively. Coded to be symmetric (could improve that easily)
    
    args of relevence:
        args.mol    - the input file (log, xyz etc) file which we will read geometry from and use to define our grid limits.
    returns:
        Tuple containing displacements suitable for the parse_dimension function called in bimolpes.py
    '''
    print(f'Estimating displacements based on final geometry of {args.mol}')
    geometry = get_geometries(args.mol)[-1]
    
    disp_xyz = split_to_xyz(args.est_disp) # get additional displacements in x ([0]), y([1]) and z([2]) from the args.est_disp argument. Parse with split_to_xyz.
    
    min_x = min_y = min_z = float('inf')    # Initialize minimums to infinity
    max_x = max_y = max_z = -float('inf')   # Initialize maxima to negative infinity
    
    for line in geometry:
        parts = line.split()
        atom_vdw_radii = Chem.GetPeriodicTable().GetRvdw(parts[0])  # get VdW radii of the atom
        x, y, z = map(float, parts[1:4])
        
        min_x = np.floor(min(min_x, 2 * (x - (atom_vdw_radii + disp_xyz[0])))) # take displacements as 2x +/- vwd radii and additional user requested displacemenet (dx/y/z)
        min_y = np.floor(min(min_y, 2 * (y - (atom_vdw_radii + disp_xyz[1]))))
        min_z = np.floor(min(min_z, 2 * (z - (atom_vdw_radii + disp_xyz[2]))))
        max_x = np.ceil(max(max_x,  2 * (x + (atom_vdw_radii + disp_xyz[0]))))
        max_y = np.ceil(max(max_y,  2 * (y + (atom_vdw_radii + disp_xyz[1]))))
        max_z = np.ceil(max(max_z,  2 * (z + (atom_vdw_radii + disp_xyz[2]))))
        
    if not args.est_hemi:
        min_z = 0
        print(f'Confining grid to z+ hemisphere (to do both +z and -z hemispheres, pass the -est_hemi flag to bimolpes.py write!')
        
    return ((min_x, max_x), (min_y, max_y), (min_z, max_z))
    
def generate_coords(extracted_geometry):
    '''
    takes extracted geometry from get_geometries and polishes it so we can reuse it in a new input .gjf file
    
    Args:
        extracted_geometry  - our extracted geometry (just a single geometry!) from get_geometries.
    
    Returns
        coords              - .gjf formatted coordinates.
    '''
    
    lines = extracted_geometry.strip().split('\n')
    data_lines = lines[5:-1]    # Skipping header and footer lines

    atomic_data = []
    coords = []                 # empty array for coordinate data

    for line in data_lines:     # Extract atomic number, element symbol, and coordinates
        columns = line.split()
        atomic_number = int(columns[1])
        element_symbol = Chem.GetPeriodicTable().GetElementSymbol(atomic_number)
        x, y, z = float(columns[3]), float(columns[4]), float(columns[5])
        atomic_data.append((element_symbol, x, y, z))

    for symbol, x, y, z in atomic_data:    # store the extracted data
        coords.append(f"{symbol}         {x} {y} {z}")
        
    return(coords)
	
def gen_grid(args):
    """
    Make a grid that spans =dx to +dx in steps of res
    
    Args:
        args    - our passed arguments
        
    args of Note:
        x/y/z   - the min/max x/y/z displacements we'll consider
        res     - the resolution (either 1, 2, or 3 values) which we'll parse with split_to_xyz
        min_dist- the minimum distance we'll allow between atoms in adjacent molecules.
        
    Returns:
        grid    - grid of points used as translation coordinates
    """
    
    res_xyz = split_to_xyz(args.res) # split our resolution up as needed
    
    x_range = np.arange(args.x[0], args.x[1] + res_xyz[0], res_xyz[0]) # here we just make a linear array of separations in x (then y, then z)
    y_range = np.arange(args.y[0], args.y[1] + res_xyz[1], res_xyz[1])
    z_range = np.arange(args.z[0], args.z[1] + res_xyz[2], res_xyz[2])

    print(f'\nGrid Limits: x={args.x} y={args.y} z={args.z}')
    print(f'Resolution:  x={res_xyz[0]} y={res_xyz[1]} z={res_xyz[2]}')
    print(f'Unpruned grid has a total of {np.prod([len(x_range), len(y_range), len(z_range)])} points')
    
    grid = [] # store the grid here
    count = 0 # count each time we don't add to the grid 
    
    for x_ in x_range: # loop over x (and y, and z) and if it passes the >= min_dist check, add this set of coords to the grid
        for y_ in y_range:
            for z_ in z_range:
                if (args.min_dist is None) or (np.sqrt(x_**2 + y_**2 + z_**2) >= args.min_dist): # min_dist check
                    grid.append((x_, y_, z_))
                else:
                    count += 1
               
    print(f'Rejected {count} grid points due to close contacts (min_dist = {args.min_dist})')
    print(f'Initial grid size: {np.shape(grid)[0]}\n')
    
    return grid

def add_fragment_label(args, geometry, fragment_number=1, atom_count = 0):
    """
    Adds "(Fragment={fragment_number})" after the atomic symbol in each line of the geometry.
    
    Args:
        args: Here we mainly want the freeze atoms 
        geometry: List of strings representing the atoms and their coordinates.
        fragment_number: the fragment number to add (e.g. 1, 2 etc.)
        atom_count: We'll use this to keep track of atom numbers between fragments 1 and 2 (because #1 of 2 is numbered as the last one in 1 + 1...)
    Returns:
        Updated geometry with "(Fragment=1)" labels as a list of strings.
        
    to do:
        need to make sure this is skipped when doing ORCA type files.
    """
    updated_geometry = []

    for x, line in enumerate(geometry):
        atom_count += 1 # add
        parts = line.split()
        if len(parts) == 4:
            atom, x, y, z = parts[0], parts[1], parts[2], parts[3]
            atom_with_fragment = f"{atom}(Fragment={fragment_number})"
            if atom_count in args.frz:
                atom_with_fragment = atom_with_fragment + ' -1  '
            if not atom_count in args.frz:
                atom_with_fragment = atom_with_fragment + ' 0   '
            new_line = f"{atom_with_fragment}   {x} {y} {z}"
            updated_geometry.append(new_line)

    return updated_geometry, atom_count


def translate_coordinates(geometry, dx=0, dy=0, dz=0):
    """
    Adjusts the geometry by adding/subtracting deltas to the x, y, and z coordinates.
    
    Args:
        geometry: List of strings representing the atoms and their coordinates.
        dx: Value to add/subtract from the x coordinate.
        dy: Value to add/subtract from the y coordinate.
        dz: Value to add/subtract from the z coordinate.
    
    Returns:
        adjusted_geometry - Adjusted geometry as a list of strings.
    """
    adjusted_geometry = []

    for line in geometry:
        parts = line.split()
        if len(parts) == 4:
            atom, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
            x_new = x + dx
            y_new = y + dy
            z_new = z + dz
            new_line = f"{atom} {x_new:.6f} {y_new:.6f} {z_new:.6f}"
            adjusted_geometry.append(new_line)

    return adjusted_geometry

def rotate_coordinates(args, coordinates, x_angle=0, y_angle=0, z_angle=0):
    '''
    Rotate coordinates by applying a rotation matrix. Takes raw geometry from get_geometry
    Args:
        coordinates - geometry from get_geometry function
        x_angle     - x_angle to rotate by, in degrees
        y_angle     - y_angle to rotate by, in degrees
        z_angle     - z_angle to rotate by, in degrees
        
    Returns:
        adjusted_geometry - the input geometry rotated by user supplied angle(s).
    '''

    rotation = R.from_euler('xyz', [x_angle, y_angle, z_angle], degrees=True) # rot matrix from SciPy Euler for brevity
    adjusted_geometry = []

    for line in coordinates:
        parts = line.split()
        if len(parts) == 4:
            atom, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
            original_coords = np.array([x, y, z])
            new_coords = rotation.apply(original_coords)
            new_line = f"{atom}    {new_coords[0]:.6f} {new_coords[1]:.6f} {new_coords[2]:.6f}"
            adjusted_geometry.append(new_line)

    return adjusted_geometry


#def distance(point1, point2):
#    """
#    Calculate the distance between two points in xyz.
#
#   This is now not used - the vectorised code should be much faster.
#    """
#    return ((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2 + (point1[2] - point2[2]) ** 2) ** 0.5
def extract_coords(fragment):
    """
    Extracts only the X, Y, Z coordinates from a molecular fragment.
    Works directly with text-based fragments (no NumPy conversions).
    """
    coords = []
    for atom in fragment:
        parts = atom.strip().split()
        if len(parts) >= 4:
            try:
                x, y, z = float(parts[-3]), float(parts[-2]), float(parts[-1])
                coords.append([x, y, z])
            except ValueError:
                raise ValueError(f"Could not extract XYZ coordinates from: {atom}")
    return np.array(coords)

def check_fragments_too_close(frag1, frag2, min_cutoff=2):
    """
    Vectorized check if any atoms in two fragments are too close.
    Uses text-based input without unnecessary pre-conversion.
    """
    coords1 = extract_coords(frag1)
    coords2 = extract_coords(frag2)

    if coords1.shape[1] != 3 or coords2.shape[1] != 3:
        raise ValueError(f"Incorrect coordinate format: frag1 {coords1.shape}, frag2 {coords2.shape}")

    distances = np.linalg.norm(coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :], axis=2)
    return np.any(distances < min_cutoff)

def check_fragments_too_far(frag1, frag2, max_cutoff=5):
    """
    Vectorized check if all atoms in two fragments are too far apart.
    Uses text-based input without unnecessary pre-conversion.
    """
    coords1 = extract_coords(frag1)
    coords2 = extract_coords(frag2)

    if coords1.shape[1] != 3 or coords2.shape[1] != 3:
        raise ValueError(f"Incorrect coordinate format: frag1 {coords1.shape}, frag2 {coords2.shape}")

    distances = np.linalg.norm(coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :], axis=2)
    return not np.any(distances <= max_cutoff)
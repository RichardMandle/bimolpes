# functions for geometry operations (geoops)

import numpy as np
from scipy.spatial.transform import Rotation as R
from rdkit import Chem

import gauops as gps

def estimate_grid_from_glog(filename, dx = 0, dy = 0, dz = 0):
    '''
    Reads a gaussian output file (glog; filename) and estimate a reasonable set of x/y/z displacements as the floor/ceiling of the min/max x/y/z values
    
    Args:
        filename: the gaussian output file from which we'll read and parse geometry
        dx/dy/dz: additional user requested displacement in X/Y/Z, respectively. Coded to be symmetric (could improve that easily)
    
    Returns:
        Tuple containing displacements suitable for the parse_dimension function called in bimolpes.py
    '''
    geometry = generate_coords(gps.get_geometries(filename)[-1])
    
    min_x = min_y = min_z = float('inf')    # Initialize minimums to infinity
    max_x = max_y = max_z = -float('inf')   # Initialize maxima to negative infinity
    
    for line in geometry:
        parts = line.split()
        atom_vdw_radii = Chem.GetPeriodicTable().GetRvdw(parts[0])  # get VdW radii of the atom
        x, y, z = map(float, parts[1:4])
        
        min_x = np.floor(min(min_x, x - (atom_vdw_radii + dx))) # take displacements +/- vwd radii and additional user requested displacemenet (dx/y/z)
        min_y = np.floor(min(min_y, y - (atom_vdw_radii + dy)))
        min_z = np.floor(min(min_z, z - (atom_vdw_radii + dz)))
        max_x = np.ceil(max(max_x, x + (atom_vdw_radii + dx)))
        max_y = np.ceil(max(max_y, y + (atom_vdw_radii + dy)))
        max_z = np.ceil(max(max_z, z + (atom_vdw_radii + dz)))
    
    return ((min_x, max_x), (min_y, max_y), (min_z, max_z))
    
def generate_coords(extracted_geometry):
    '''
    takes extracted geometry from get_geometries and polishes it so we can reuse it in a new input .gjf file
    
    Args:
    extracted_geometry - our extracted geometry (just a single geometry!) from get_geometries.
    
    Returns
    coords - .gjf formatted coordinates.
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
	
def gen_grid(x, y, z, res = 0.25):
    """
    Make a grid that spans =dx to +dx in steps of res
    
    Args:
        x,y,z   - tuple; min/max x-value (and y, z)
        res     - resolution
    Returns:
        grid    - grid of points used as translation coordinates
    """
        
    x_range = np.arange(x[0], x[1] + res, res)
    y_range = np.arange(y[0], y[1] + res, res)
    z_range = np.arange(z[0], z[1] + res, res)

    grid = []
    for x_ in x_range:
        for y_ in y_range:
            for z_ in z_range:
                grid.append((x_, y_, z_))
    print(f'Using a grid of x={x} y={y} z={z} res={res}; Final grid size={np.product(np.shape(grid))}')
    
    return grid

def add_fragment_label(geometry, fragment_number=1):
    """
    Adds "(Fragment={fragment_number})" after the atomic symbol in each line of the geometry.
    
    Args:
        geometry: List of strings representing the atoms and their coordinates.
    
    Returns:
        Updated geometry with "(Fragment=1)" labels as a list of strings.
    """
    updated_geometry = []

    for line in geometry:
        parts = line.split()
        if len(parts) == 4:
            atom, x, y, z = parts[0], parts[1], parts[2], parts[3]
            atom_with_fragment = f"{atom}(Fragment={fragment_number})"
            new_line = f"{atom_with_fragment} {x} {y} {z}"
            updated_geometry.append(new_line)

    return updated_geometry


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


def rotate_coordinates(coordinates, x_angle=0, y_angle=0, z_angle=0):
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
    x_rad = np.radians(x_angle)
    y_rad = np.radians(y_angle)
    z_rad = np.radians(z_angle)

    rotation = R.from_euler('xyz', [x_angle, y_angle, z_angle], degrees=True) # rot matrix from SciPy Euler for brevity
    adjusted_geometry = []

    for line in coordinates:
        parts = line.split()
        if len(parts) == 4:
            atom, x, y, z = parts[0], float(parts[1]), float(parts[2]), float(parts[3])
            original_coords = np.array([x, y, z])
            new_coords = np.dot(R, original_coords)
            new_line = f"{atom}    {new_coords[0]:.6f} {new_coords[1]:.6f} {new_coords[2]:.6f}"
            adjusted_geometry.append(new_line)

    return adjusted_geometry

def distance(point1, point2):
    """
    Calculate the distance between two points in xyz.
    """
    return ((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2 + (point1[2] - point2[2]) ** 2) ** 0.5

def check_fragments_too_close(frag1, frag2, min_cutoff=2):
    """
    Check if any atoms in two fragments are too close based on the cutoff distance.

    Args:
        frag1 (list): The list of atoms and their coordinates for the first fragment.
        frag2 (list): The list of atoms and their coordinates for the second fragment.
        cutoff (float): The cutoff distance to consider atoms too close.

    Returns:
        bool: True if any atoms are too close, False otherwise.
    """
    for atom1 in frag1:
        for atom2 in frag2:
            coords1 = tuple(map(float, atom1.split()[1:4]))
            coords2 = tuple(map(float, atom2.split()[1:4]))
            if distance(coords1, coords2) < min_cutoff:
                return True
    return False

def check_fragments_too_far(frag1, frag2, max_cutoff=5):
    """
    Check if any atoms in two fragments are too close based on the cutoff distance.

    Args:
        frag1 (list): The list of atoms and their coordinates for the first fragment.
        frag2 (list): The list of atoms and their coordinates for the second fragment.
        cutoff (float): The cutoff distance to consider atoms too close.

    Returns:
        bool: True if any atoms are too close, False otherwise.
    """
    for atom1 in frag1:
        for atom2 in frag2:
            coords1 = tuple(map(float, atom1.split()[1:4]))
            coords2 = tuple(map(float, atom2.split()[1:4]))
            if distance(coords1, coords2) <= max_cutoff:
                return False 
    return True  

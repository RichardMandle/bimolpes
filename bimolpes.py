# BIMOLPES
# Dr. R.J.Mandle; UoL, 2024.

import numpy as np
import platform
import argparse
import re
import os
import glob

# import our own modules
import gauops as gps
import geoops as geo
import visualisation as vis
import processing as pro

def main():
    parser = argparse.ArgumentParser(description="bimolpes - generate and process bimolecular potential energy surfaces")
    subparsers = parser.add_subparsers(help='commands', dest='command') # subparsers for write and read modes
    subparsers.required = True

    # Sub-parser for the write command
    write_parser = subparsers.add_parser('write', help='Reads a Gaussian output file (or 2 files) and extracts molecular geometry\nGenerates and writes a grid of translated coordinates for two molecules as .gjf files, including job sub script (.sh for SGE)\n\n Files are given as .zip; just upload to HPC unzip and run')
    write_parser.add_argument("-x", type=parse_dimension, default=(0.0, 0.0), help="Range for x dimension, specify as 'start end' or just 'value' to use '-value value'")
    write_parser.add_argument("-y", type=parse_dimension, default=(0.0, 0.0), help="Range for y dimension, specify as 'start end' or just 'value' to use '-value value'")
    write_parser.add_argument("-z", type=parse_dimension, default=(0.0, 0.0), help="Range for z dimension, specify as 'start end' or just 'value' to use '-value value'")
    write_parser.add_argument('-est', action='store_const', const=True, default=True, help='Estimate grid dimensions from input file (experimental)')
    write_parser.add_argument("-est_x", type=float, default=0.0, help="additional displacement in x dimension from estimated displacement")
    write_parser.add_argument("-est_y", type=float, default=0.0, help="additional displacement in y dimension from estimated displacement")
    write_parser.add_argument("-est_z", type=float, default=0.0, help="additional displacement in z dimension from estimated displacement")
    write_parser.add_argument('-res', type=float, default=1.0, help='Resolution of the grid')    
    
    write_parser.add_argument('-xa', type=float, default=0.0, help='Angle to rotate mol. #2 by in X-axis, in degrees')
    write_parser.add_argument('-ya', type=float, default=0.0, help='Angle to rotate mol. #2 by in Y-axis, in degrees')
    write_parser.add_argument('-za', type=float, default=0.0, help='Angle to rotate mol. #2 by in Z-axis, in degrees')
    
    write_parser.add_argument('-min_dist', type=float, default=3.0, help='Minimum intermolecular separation of molecules in grid, in Å')
    write_parser.add_argument('-max_dist', type=float, default=6.0, help='Maximum intermolecular separation of molecules in grid, in Å')
    
    write_parser.add_argument('-mem', type=int, default=4, help='RAM for gaussian job, in GB (default = 4)')    
    write_parser.add_argument('-cpu', type=int, default=4, help='CPU cores for gaussian job (default = 4)')     
    write_parser.add_argument('-groute', type=str, default='B3LYP cc-pVTZ EmpiricalDispersion=GD3BJ counterpoise=2', help='Specify the full Gaussian route section; default B3LYP cc-pVTZ EmpiricalDispersion=GD3BJ Counterpoise=2')     
        
    write_parser.add_argument('-inp', type=str, default='', help='Gaussian .log file of OPTIMISED structure of molecule 1')    
    write_parser.add_argument('-inp2', type=str, default=None, help='Gaussian .log file of OPTIMISED structure of molecule 2 (optional, defaults to homo-bimolecular PES')    
    write_parser.add_argument('-out', type=str, default=None, help='Naming schema for output .gjf files (optional; defaults to mol1_mol2\mol1_mol2 <<TODO extending iteratively?>>)')   
    write_parser.add_argument('-zip', action='store_const', const=False, default=True, help='Don\'t zip output!')
    
    write_parser.set_defaults(func=write_grid)

    # Sub-parser for the read command; first up, options for the 3d scatter plot.
    read_parser = subparsers.add_parser('read', help='Read completed gaussian calculations generated using write mode\nVisualise data, molecular contacts, finds minima')
    read_parser.add_argument('-path', type=str, default=os.getcwd(), help='Path to read .log files from')
    read_parser.add_argument('-method', type=int, default=1, help='Energy type to use; 0 = SCF, 1 = counterpoise corrected complexation')

    read_parser.add_argument('-nosave', action='store_const', const=True, default=False, help='Save output to .npz for easy reloading')
    read_parser.add_argument('-reload', action='store_const', const=True, default=False, help='reload bimolpes data from .npz file')
    read_parser.add_argument('-filename', type=str, default='mydata', help='filename to use for saving .npz')  

    read_parser.add_argument('-minima', type=int, default=10, help='Number of minima sites to consider')
    read_parser.add_argument('-ethr', type=int, default=8, help='Energy cutoff for identifying discrete minima')
    read_parser.add_argument('-dthr', type=int, default=1, help='Distance cutoff for identifying discrete minima')           
    read_parser.set_defaults(func=handle_read)
    
    # Sub-parser for the plot command; first up, options for the 3d plot.
    plot_parser = subparsers.add_parser('plot', help='Plotting and visualising data using mayavi')
    plot_parser.add_argument('-filename', type=str, default='mydata', help='filename to use for saving .npz')  
    
    plot_parser.add_argument('-plt_emax', type=float, default=0, help='Only show points where energy is below e_max (float)')    
    plot_parser.add_argument('-plt_size', type=float, default=1, help='Size of grid energy points (float, default = 1)')    
    plot_parser.add_argument('-plt_alpha', type=float, default=1, help='Alpha vlaues of grid energy points (float, default = 1)')    
    plot_parser.add_argument('-plt_style', type=str, default='surface', help='Select between \'surface\', \'points\', \'wireframe\', \'fancymesh\' etc; see mlab')
    plot_parser.add_argument('-plt_line', type=float, default=2.0, help='float; sets the line_width for some types of plot, e.g. \'wireframe\', \'fancymesh\'; see mlab docs')
    
    plot_parser.add_argument('-plt_cmap', type=str, default='plasma', help='Pass any valid matplotlib colourmap (cmap) to use in plotting')
                
    # sub parts for drawing of molecule using matplotlib...
    plot_parser.add_argument('-mol', type=str, default=None, help='Gaussian .log file for the molecule to draw on the PES')
    plot_parser.add_argument('-mol2', type=str, default=None, help='Gaussian .log file for the molecule to draw on the PES')
    plot_parser.add_argument('-mol2_idx', type=int, default=0, help='Draw a second molecule on the PES with the given index value')   
    plot_parser.add_argument('-real_size', action='store_const', const=False, default=True, help='use more "realistic" atom sizes when drawing')
    plot_parser.add_argument('-mol_size', type=float, default=2, help='Size of \'atoms\' in molecule (float, default = 128)')    
    plot_parser.add_argument('-mol_alpha', type=float, default=1, help='Alpha vlaues of atoms in molecule (float, default = 1)')    
    plot_parser.add_argument('-greyscale', action='store_const', const=True, default=False, help='Draw the molecule(s) in greyscale')

    plot_parser.set_defaults(func=plot_data)

    args = parser.parse_args()

    args.func(args)
    
def plot_data(args):
    # master function for plotting using mayavi
    from mayavi import mlab
    data = pro.reload_data_from_npz(args.filename)
    dx_values, dy_values, dz_values, dyz_values, e_values, de_values = pro.shape_data(data)
    
    if args.plt_emax != 0:
        dx_values, dy_values, dz_values, e_values = pro.mask_values(dx_values, dy_values, dz_values, e_values, args.plt_emax)
        
    fig = mlab.figure(bgcolor=(1, 1, 1), size=(800, 600))
    vis.plot_volumetric_mayavi(dx_values, dy_values, dz_values, e_values, e_max=args.plt_emax, size=args.plt_size,  line_width = args.plt_line, alpha=args.plt_alpha, plot_style=args.plt_style, cmap=args.plt_cmap, fig = fig)

    if args.mol != None:
        vis.plot_molecule_mayavi(filename=args.mol, size = args.mol_size, alpha = args.mol_alpha, real_size = args.real_size, greyscale = args.greyscale, fig = fig)

    if args.mol2_idx != 0:
        if args.mol2 == None:
            filename = args.mol
        if args.mol2 != None:
            filename = args.mol2
            
        vis.plot_molecule_mayavi(filename=args.mol2, displacements = data[args.mol2_idx,0:3], size = args.mol_size, alpha = args.mol_alpha, real_size = args.real_size, greyscale = args.greyscale, fig = fig)
    mlab.show() 
        
# functions for generating grid:
def parse_dimension(arg_value):
    '''
    Parsing function; splits input allowing many delimiters, returns x/y/z limits
    '''
    values = re.split(r'[^\d.-]+', arg_value)
    values = list(map(float, filter(None, values)))
    if len(values) == 1:
        return (- values[0], values[0])  # Default to range from value to value if one is given
    elif len(values) == 2:
        return tuple(values)
    else:
        raise argparse.ArgumentTypeError(f"Expected 1 or 2 values, got {len(values)} in input {arg_value}")

def write_grid(args):
    if args.out == None:
        if args.inp2 == None:
            outpath = (args.inp).split('.')[0]
        if (args.inp2 != None) & (args.inp != None):
            outpath = (args.inp).split('.')[0] + '_' + (args.inp2).split('.')[0]
        
    if args.out != None:
        outpath = args.out
        
    os.makedirs(outpath, exist_ok=True)
    file_plus_path = outpath + '\\' + outpath
        
    if args.est == True:
        print(f'Estimating displacements based on final geometry of {args.inp}')
        args.x, args.y, args.z = geo.estimate_grid_from_glog(args.inp, dx = args.est_x, dy = args.est_y, dz = args.est_z)          # here, estimate the x/y/z ranges from the input file using geo.estimate_grid_from_glog.
        
    grid = geo.gen_grid(x = args.x, y = args.y, z = args.z, res = args.res, min_disp = args.min_dist)     # get a grid using the user supplied spacings.
    
    outname, count =  create_grid(args.inp, grid = grid, x_angle = args.xa, y_angle = args.ya, z_angle = args.za, min_dist = args.min_dist, max_dist = args.max_dist, filename2 = args.inp2, outname = file_plus_path)
    
    gps.make_sge_job(filename = file_plus_path, nproc = args.cpu, vmem = args.mem, startjob=1, endjob=count) # make the SGE job

    if args.zip == True: # bundle into .zip for unleashing on HPC (quicker to upload one file than 100k small ones)
        pro.zip_files(dir_path = outpath, zip_file = outpath, ext =('gjf','sh'))
    
def create_grid(filename, grid, nproc=8, vmem=4, x_angle = 0, y_angle = 0, z_angle = 0, min_dist = 2, max_dist = 8, count = 1, filename2 = None, outname = None):                
    '''                                                                                                                   
    Function to create a grid of points which we use as translation vectors for our second molecule                       
    For each point on the grid we offset molecule 2 by this much, allowing us to scan across all 3 dimensions             
                                                                                                                          
    Args:                                                                                                                 
        grid      - grid generated by geo.gen_grid()
        
        rotate    - Bool; do rotation or not? Defaults to True if angle is supplied
        x_angle   - rotation about x-axis (degrees)
        y_angle   - rotation about y-axis (degrees)
        z_angle   - rotation about z-axis (degrees)

        min_dist  - cutoff if any atom on mol1 is this close to any atom in mol2 then this set of coordinates is rejected
        max_dist  - cutoff if any atom on mol1 is this far to any atom in mol2 then this set of coordinates is rejected
    
        count     - a counter used for writing the right number of SGE files (default = 1)
    
        filename2 - Can supply a different molecule for the second one to "scan"
    Returns:
        gjf files - Gaussian input files, one for each accepted set of translations
        sh file   - SGE file for executing on ARC3/4
        
    '''
    
    outname = pro.get_output_name(outname, filename, filename2)
       
    gps.clean_path(outname) # tidy directory
    
    geometry = geo.generate_coords(gps.get_geometries(filename)[-1])
    frag2_geometry = geometry    # set as the same for now, update below

    if (filename2 != None) | ((filename2 != filename) & (filename2 != None)):
        frag2_geometry = geo.generate_coords(gps.get_geometries(filename2)[-1])
        
    if (np.abs(x_angle) + np.abs(y_angle) + np.abs(z_angle)) != 0:
        frag2_geometry = geo.rotate_coordinates(geometry,  x_angle=x_angle, y_angle=y_angle, z_angle=z_angle)

    frag1 = geo.add_fragment_label(geometry,1)
    
    for gri in grid:
        frag2 = geo.add_fragment_label(geo.translate_coordinates(frag2_geometry, dx=gri[0], dy=gri[1], dz=gri[2]), 2)

        if not geo.check_fragments_too_close(frag1, frag2, min_cutoff=min_dist):
            if not geo.check_fragments_too_far(frag1, frag2, max_cutoff=max_dist):
                displacement = f'dx={gri[0]}/dy={gri[1]}/dz={gri[2]}/'
                file_name = outname + '_' + str(count)
                gps.write_gjf(frag1, frag2, displacement=displacement, file_name=file_name)
                count += 1
    print(f'A total of {count-1} geometries written to .gjf after spatial cutoffs of min_dist={min_dist} and max_dist={max_dist}')
    
    return outname, count
    

def handle_read(args):
    # master function; reads data, saves data, prints info on minima, .
   
    if args.reload == True:
        data = pro.reload_data_from_npz(args.filename)
    else:
        data = pro.process_files(path=args.path, method=args.method)
    
    if args.nosave == False:
        pro.save_data_to_npz(args.filename, data)
        
    dx_values, dy_values, dz_values, dyz_values, e_values, de_values = pro.shape_data(data)
       
    minima_list = pro.find_local_minima(data, num_minima = args.minima, e_threshold = args.ethr , d_threshold = args.dthr)
    
    minima_indices = pro.find_indices_of_minima(data, minima_list)
    
    print ("\n{:<12} {:<7} {:<7} {:<7} {:<7} {:<7} ".format('Type', 'X pos.', 'Y pos.', 'Z pos.', 'Energy', 'Index'))
    for i, minima in enumerate(minima_list):
        if i == 0:
            print ("{:<12} {:<7} {:<7} {:<7} {:<7} {:<7}".format('Global Min.', minima[0],minima[1],minima[2],minima[3],minima_indices[i]))
        if i != 0:
            print ("{:<12} {:<7} {:<7} {:<7} {:<7} {:<7}".format('Local Min.', minima[0],minima[1],minima[2],minima[3],minima_indices[i]))

if __name__ == "__main__":
    print(' ___ ___ __  __  ___  _    ___ ___ ___') 
    print('| _ )_ _|  \/  |/ _ \| |  | _ \ __/ __|')
    print('| _ \| || |\/| | (_) | |__|  _/ _|\__ \\')
    print('|___/___|_|  |_|\___/|____|_| |___|___/')
    print(f'\nVersion: 0.6; running on {platform.system()} {platform.version()}')
    print(f'Authors: Dr. R.J.Mandle; University of Leeds, 2024\n')
    main()
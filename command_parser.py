import argparse
import os
import re

class BimolPESParser:
    '''
    Parser class to handle input arguments for bimolpes.
    '''
    def __init__(self, command_functions, config = None):
        '''
        Args:
            command_functions   - passed list of command functions from bimolpes.py (write, read, plot)
            config              - configuration loaded from config.ini by bimolpes.py
        '''
        self.command_functions = command_functions
        self.config = config
        self.parser = argparse.ArgumentParser(description="bimolpes - generate and process bimolecular potential energy surfaces")
        self.subparsers = self.parser.add_subparsers(help='Specify either write, read, or plot modes.\n\n write - builds a grid of translated coordinates (e.g. calculation setup\n read  - reads DFT output and saves data to .npz for easy reloading/plotting\n plot  - plots data using mayavi', dest='command', required=True)
        self.setup_write_parser()
        self.setup_read_parser()
        self.setup_plot_parser()
        
    def setup_write_parser(self):
        write_parser = self.subparsers.add_parser('write', help='Reads a Gaussian output file (or 2 files) and extracts molecular geometry\nGenerates and writes a grid of translated coordinates for two molecules as .gjf files, including job sub script (.sh for SGE)\n\n Files are given as .zip; just upload to HPC unzip and run')
        write_parser.add_argument("-x", type=parse_dimension, default=(0.0, 0.0), help="Range for x dimension, specify as 'start end' or just 'value' to use '-value value'")
        write_parser.add_argument("-y", type=parse_dimension, default=(0.0, 0.0), help="Range for y dimension, specify as 'start end' or just 'value' to use '-value value'")
        write_parser.add_argument("-z", type=parse_dimension, default=(0.0, 0.0), help="Range for z dimension, specify as 'start end' or just 'value' to use '-value value'")
        write_parser.add_argument('-est', action='store_true', default=self.config.getboolean('GRID','EstimateGrid', fallback = False), help='Estimate grid dimensions from input file (experimental)')
        write_parser.add_argument('-est_hemi', action='store_true', default=self.config.getboolean('GRID','EstimateGridHemisphere', fallback = False), help='When estimating grid, do both hemispheres (i.e. don\'t stop at z = 0!')
        write_parser.add_argument("-est_disp", type=parse_resolution, default=self.config.get('GRID','AdditionalDisplacement', fallback = '0.0'), help="additional displacement from estimated grid;  provide one, two, or three values for x, y, and z respectively")
        write_parser.add_argument('-res', type=parse_resolution, default=self.config.get('GRID','Resolution', fallback = '1.0'), help='Resolution of the grid; provide one, two, or three values for x, y, and z respectively')
        write_parser.add_argument('-rot', type=parse_resolution, default=self.config.get('GRID','Rotation', fallback = '0.0,0.0,0.0'), help='Angle to rotate mol. #2 by in in degrees; one, two, or three values for rotation about x, y, and z respectively ')
        write_parser.add_argument('-backoff', action='store_true', default=self.config.getboolean('GRID','TurnOffBackUp', fallback = False), help='Turn off backing up data - default is to backup a given output directory to "directory.zip" rather than overwrite all the contents')
              
        write_parser.add_argument('-min_dist', type=float, default=self.config.getfloat('GRID','MinDist', fallback = 3.0), help='Minimum intermolecular separation of molecules in grid, in Å (fallback = 3.0)')
        write_parser.add_argument('-max_dist', type=float, default=self.config.getfloat('GRID','MaxDist', fallback = 6.0), help='Maximum intermolecular separation of molecules in grid, in Å (fallback = 6.0)')
        
        write_parser.add_argument('-mem', type=int, default=self.config.getint('GAUSSIAN','RAM', fallback = 8), help='RAM for gaussian job, in GB (int; fallback = 8)')    
        write_parser.add_argument('-cpu', type=int, default=self.config.getint('GAUSSIAN','CPU', fallback = 8), help='CPU cores for gaussian job (int; fallback = 8)')     
        write_parser.add_argument('-groute', type=str, default=self.config.get('GAUSSIAN','GRoute'), help='Specify the full Gaussian route section; default: #T B3LYP cc-pVTZ EmpiricalDispersion=GD3BJ Counterpoise=2')     
        write_parser.add_argument('-gver', type=str, default=self.config.get('GAUSSIAN','Version', fallback = 'G16'), help='Specify Gaussian version used for calculations (str; fallback = G16)') 
        write_parser.add_argument('-disk', type=str, default=self.config.get('GAUSSIAN','MaxDisk', fallback = '5'), help='Specify Max Disk option for Gaussian calculations in GB (str; fallback = 5)') 
               
        write_parser.add_argument('-inp', type=str, default='', help='Gaussian .log file of OPTIMISED structure of molecule 1')    
        write_parser.add_argument('-inp2', type=str, default=None, help='Gaussian .log file of OPTIMISED structure of molecule 2 (optional, defaults to homo-bimolecular PES')    
        write_parser.add_argument('-out', type=str, default=None, help='Naming schema for output .gjf files (optional; defaults to mol1_mol2\mol1_mol2 <<TODO extending iteratively?>>)')   
        write_parser.add_argument('-zip', action='store_const', const=False, default=True, help='Don\'t zip output!')
        
        write_parser.add_argument('-frz', type=parse_atom_numbers, default='', help='Int; Specify atom numbers to freeze during the calculation. This is useful if also optimising; pass a set of 3 or more non-colinear points to freeze a molecule but allow (some) internal bond rotations etc. Can be done on molecule #1 and #2; be mindful that the atom numbers need to be set by checking the number of atoms in the input .log file(s). Allows partial optimisation (a semi-relaxed scan?) i.e. molecules frozen in place centre-to-centre but most atoms/bonds allowed to relax.')
        write_parser.set_defaults(func=self.command_functions['write_grid'])

    def setup_read_parser(self):
        # Sub-parser for the read command; first up, options for the 3d scatter plot.
        read_parser = self.subparsers.add_parser('read', help='Read completed gaussian calculations generated using write mode\nVisualise data, molecular contacts, finds minima')
        read_parser.add_argument('-path', type=str, default=os.getcwd(), help='Path to read .log files from')
        read_parser.add_argument('-method', type=int, default=self.config.getint('READ','Method', fallback = 1), help='Energy type to use; 0 = SCF, 1 = counterpoise corrected complexation')

        read_parser.add_argument('-nosave', action='store_const', const=False, default=False, help='Don\t save output to .npz')
        read_parser.add_argument('-reload', action='store_const', const=True, default=False, help='reload bimolpes data from .npz file')
        read_parser.add_argument('-filename', type=str, default='mydata', help='specify filename to use for saving .npz data; (default = mydata)')  

        read_parser.add_argument('-minima', type=int, default=self.config.getint('READ','Minima', fallback = 10), help='Number of minima sites to consider (int; fallback = 1)')
        read_parser.add_argument('-ethr', type=float, default=self.config.getfloat('READ','EThr', fallback = 8.0), help='Energy cutoff for identifying discrete minima (float; fallback = 8.0')
        read_parser.add_argument('-dthr', type=float, default=self.config.getfloat('READ','DThr', fallback = 1.0), help='Distance cutoff for identifying discrete minima (float; fallback = 1.0)')           
        
        read_parser.set_defaults(func=self.command_functions['handle_read'])
        
    def setup_plot_parser(self):
        # Sub-parser for the plot command; first up, options for the 3d plot.
        plot_parser = self.subparsers.add_parser('plot', help='Plotting and visualising data using mayavi')
        plot_parser.add_argument('-filename', type=str, default='mydata', help='filename to use for saving .npz')  
        
        plot_parser.add_argument('-plt_emax', type=float, default=self.config.getfloat('PLOT','PltEMax', fallback = 0), help='Only show points where energy is below e_max (float; fallback = 0)')    
        plot_parser.add_argument('-plt_size', type=float, default=self.config.getfloat('PLOT','PltSize', fallback = 1), help='Size of grid energy points (float; default = 1)')    
        plot_parser.add_argument('-plt_alpha', type=float, default=self.config.getfloat('PLOT','PltAlpha', fallback = 0.25), help='Alpha vlaues of grid energy points (float; default = 1)')    
        plot_parser.add_argument('-plt_style', type=str, default=self.config.get('PLOT', 'PltStyle', fallback='points'), help='Select between \'surface\', \'points\', \'wireframe\', \'fancymesh\' etc; see mlab (str; fallback = points)')
        plot_parser.add_argument('-plt_line', type=float, default=self.config.getfloat('PLOT','LineWidth', fallback = 2.0), help='Sets the line_width for some types of plot, e.g. \'wireframe\', \'fancymesh\'; see mlab docs (float; fallback = 2.0)')
        plot_parser.add_argument('-plt_cmap', type=str, default=self.config.get('PLOT', 'ColorMap', fallback='viridis'), help='Pass any valid matplotlib colourmap (cmap) to use in plotting (str; fallback = viridis)')
        plot_parser.add_argument('-plt_res', type=float, default=self.config.get('PLOT', 'Resolution', fallback=20.0), help='Resolution of points/surface to use in plotting (float; fallback = 20.0)')
        plot_parser.add_argument('-plt_flipz', action='store_true', default=self.config.getboolean('PLOT','FlipZ', fallback = False), help='Flip data about the Z-axis; effectively places the molecule "on top" of the surface data, just a visual effect for publications.')
        plot_parser.add_argument('-plt_max', type=float, default=self.config.get('PLOT', 'VMax', fallback=0), help='VMax of points/surface to use in plotting (float; fallback = 0)')
        plot_parser.add_argument('-plt_min', type=float, default=self.config.get('PLOT', 'VMin', fallback=0), help='VMin of points/surface to use in plotting (float; fallback = 0)')    
        
        # sub parts for drawing of molecule
        plot_parser.add_argument('-mol', type=str, default=None, help='Gaussian .log file for the molecule to draw on the PES')
        plot_parser.add_argument('-mol2', type=str, default=None, help='Gaussian .log file for the molecule to draw on the PES')
        plot_parser.add_argument('-mol2_idx', type=int, default=0, help='Draw a second molecule on the PES with the given index value')   
        plot_parser.add_argument('-mol_size', type=float, default=self.config.getfloat('PLOT','MolSize', fallback = 2), help='Size of \'atoms\' in molecule (float, default = 128)')    
        plot_parser.add_argument('-mol_real_size', action='store_true', default=self.config.getboolean('PLOT','RealSize', fallback = False),  help='use more "realistic" atom sizes when drawing')
        plot_parser.add_argument('-mol_alpha', type=float, default=self.config.getfloat('PLOT','MolAlpha', fallback = 1), help='Alpha vlaues of atoms in molecule (float, default = 1)')    
        plot_parser.add_argument('-mol_grey', action='store_true', default=self.config.getboolean('PLOT','Greyscale', fallback = False), help='Draw the molecule(s) in greyscale')
        plot_parser.add_argument('-mol_white', action='store_true', default=self.config.getboolean('PLOT','ColorWhite', fallback = False), help='Draw the molecule(s) with all atoms coloured white')
        
        # sub parts used for making fairly specific figures
        plot_parser.add_argument('-vis_plt', action='store_const', const=True, default=False, help='if true, trigger doing a vis_plt')
        plot_parser.add_argument("-vis_plt_x", type=parse_dimension, default=(0.0, 0.0), help="Range for x dimension, specify as 'start end' or just 'value' to use '-value value'")
        plot_parser.add_argument("-vis_plt_y", type=parse_dimension, default=(0.0, 0.0), help="Range for y dimension, specify as 'start end' or just 'value' to use '-value value'")
        plot_parser.add_argument("-vis_plt_z", type=parse_dimension, default=(0.0, 0.0), help="Range for z dimension, specify as 'start end' or just 'value' to use '-value value'")
        plot_parser.add_argument("-vis_plt_size", type=float, default=2.0, help="Range for z dimension, specify as 'start end' or just 'value' to use '-value value'")

        plot_parser.set_defaults(func=self.command_functions['plot_data'])
        
    def get_parser(self):
        '''
        Attempt to gracefully handle errors if strange things are passed as arguments.
        '''
        try:
            return self.parser
        except argparse.ArgumentError as e:
            print(f"Argument error: {e}")
            self.parser.print_help()


            
def parse_atom_numbers(num_value):
    '''
    Parsing function that takes a string of numbers separated by various delimiters
    and splits them into a list of floats.
    
    Args:
        num_value (str): A string containing numbers separated by delimiters.
    
    Returns:
        list of float: A list of parsed numbers as floats.
    '''
    
    values = re.split(r'[^\d.-]+', num_value)
    
    return list(map(int, filter(None, values)))
            
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

def parse_resolution(res_value):
    '''
    Parse the input to allow one to three resolution values
    '''
    res_value = re.split(r'[^\d.-]+', res_value)
    res_values = list(map(float, res_value))
    if len(res_values) > 3:
        raise argparse.ArgumentTypeError("Maximum of three resolution values allowed (x, y, z).")
    return res_values
# BIMOLPES
# Dr. R.J.Mandle; UoL, 2024.

import platform
import configparser
import os

# import our own modules
import geoops as geo
import processing as pro
import visualisation as vis
    
# import the BimolPESParser class for handling input arguments.
from command_parser import BimolPESParser

def main():
    '''
    we define the command_functions here (and call them in command_parser! so that we can
    keep the business end of the code seperate from the parser class.
    
    To add new command_functions, update the dict and add new entries to the Parser class...

    Put functions in the appropriate modules! 
    '''
    
    config = configparser.ConfigParser()
    config_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),'config.ini')
    config.read(config_file)
    
    command_functions = {
        'write_grid': geo.write_grid,
        'handle_read': pro.handle_read,
        'plot_data': vis.plot_data
    }
    
    parser = BimolPESParser(command_functions, config).get_parser()
    args = parser.parse_args()

    if hasattr(args, 'func'):
        args.func(args)
    else: # help them out if needed
        parser.print_help()
                
if __name__ == "__main__":
    print(' ___ ___ __  __  ___  _    ___ ___ ___') 
    print('| _ )_ _|  \/  |/ _ \| |  | _ \ __/ __|')
    print('| _ \| || |\/| | (_) | |__|  _/ _|\__ \\')
    print('|___/___|_|  |_|\___/|____|_| |___|___/')
    print(f'\nVersion: 0.8; running on {platform.system()} {platform.version()}')
    print(f'Authors: Dr. R.J.Mandle; University of Leeds, 2024\n')
    main()
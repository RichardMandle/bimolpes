# functions for gaussian operations (gauops)

import numpy as np
import re
import os
import glob

import geoops as geo

def get_geometries(file_path):
    '''
    basic function that extracts any and all geometries present in a Gaussian output file.
    
    These can be the output of an optimisation (e.g. all_geometries[-0], the last geometry)
    They might be all geometries in an optimisation of some sort
    
    Args:
        file_path       - the file you want to read
    
    Returns:
        all_geometries  - the extracted geometry data from the "standard orientation" tables.
    '''

    all_geometries = [] # List to store the extracted lines

    with open(file_path, 'r') as file:
        file_text = file.read()
        match_start = [m.start() for m in re.finditer('Standard orientation',file_text)]
        match_finish = [n.start() for n in re.finditer('Rotational constants',file_text)]

    for m in range(len(match_start)):
        all_geometries.append(file_text[match_start[m]:match_finish[m]])

    return all_geometries

def read_gaussian_route_section(file_path):
    '''
    Reads the route-section of a Gaussian .log/.out file
    and returns the route options as a list
    
    Args:
        file_path - file path of the .log/.out file to read
        
    Returns:
        input_options - a list of input options
    '''
    
    input_options=[]
    with open(file_path, 'r') as log_file:
        start_flag = False
        for line in log_file:
            line = line.strip()
            if line.startswith('----------------------------------------------------------------------'):
                if start_flag:
                    break  # Stop processing after the first set of input options
                start_flag = True
            elif start_flag:
                options = line.split()
                input_options.extend(options)

    return(input_options) 
	
def write_gjf(args, 
              frag1, 
              frag2,
              displacement,
              file_name='test'):
    '''
    Tool for writing two molecular geometries (frag1, frag2) into a gaussian .gjf input file.
    
    Args:
        args            - arguments from bimolpes.py; contains info on CPU cores, memory, maxdisk, Gaussian route.
        frag1(2)        - fragment 1 (2) XYZ information formatted as per .gjf file with fragment=XYZ
        displacement    - job name, but with displacement coordinates encoded as dx/y/z and rotation as rx/y/z
        file_name       - output file name of .gjf file
        
    Discussion:
        The spin/charge parameter is hardcoded as '0 1'; should this be '0 1 0 1 0 1', as we have a bimolecular .gjf?
        It might also be worth allowing the user to pass a custom spin/charge here.
    '''
    
    def format_coordinate(coord):
        try:
            float_val = float(coord)
            if float_val.is_integer():
                return "{:d}".format(int(float_val))
            else:
                return "{:.6e}".format(float_val)
        except ValueError:
            return str(coord)

    with open(file_name + '.gjf', 'w') as f:
        f.write(f'%nprocshared={args.cpu}\n')
        f.write(f'%mem={args.mem}GB\n')
        
        if args.chk:  # write our .chk checkpoint if requested
            chk_filename = os.path.basename(file_name)  # get the base name of the file
            f.write(f'%chk={chk_filename}.chk\n')
            
        f.write(f'{args.groute} maxdisk={args.disk}GB\n\n') 
        f.write(displacement + '\n\n')
        f.write('0 1\n') # Here we provide spin/charge info; might want to allow flexibility here

        for atom in frag1 + frag2:
            parts = atom.split()
            formatted_parts = [parts[0], parts[1]] + [format_coordinate(coord) for coord in parts[2:]]
            formatted_atom = " ".join(formatted_parts)
            f.write(formatted_atom)
            f.write('\n')
        f.write('\n\n')
    return file_name
    
def make_sge_job(args, outname, startjob=0, endjob=0):
    """
    Generate a job script for running a Gaussain job on ARC (the UoL compute clusters).
    
    Args:
        args:           command line arguments passed from elsewhere.
        outname (str): The name of the job script file (default: 'noname').
        startjob (int): The starting index of the task array (default: 0).
        endjob (int): The ending index of the task array (default: 0).

    Returns:
        None
    """
    file_name, _ = os.path.splitext(os.path.basename(outname))

    # Decide if we have a task array (i.e. multiple .gjf files) 
    if startjob == endjob:
        multiple = False
        
    if startjob != endjob:
        multiple = True
        
    with open(outname + '.sh', 'w') as f:
        f.write('#$ -cwd \n')
        f.write('#$ -V\n')
        f.write('#$ -l h_rt=48:00:00\n')
        f.write(f'#$ -l h_vmem={args.mem}G\n')
        f.write(f'#$ -pe smp {args.cpu}\n')
        f.write(f'#$ -l disk={args.disk}G\n')
        
        if multiple:
            f.write(f'#$ -t {startjob}-{endjob}\n')  # Create a task array

        f.write('module add gaussian\n')
        f.write('export GAUSS_SCRDIR=$TMPDIR\n')
        
        f.write(f'{args.gver} {file_name}{"_$SGE_TASK_ID" if multiple else ""}.gjf\n')
        
        f.write('rm *core* *.sh.* \n') #if job fails, this line cleans up messy core files
        
        if args.chk: # if requested checkpoint we likely want the formatted checkpoint too; so do this:
            chk_filename = os.path.basename(file_name)
            f.write(f'formchk {chk_filename}{"_$SGE_TASK_ID" if multiple else ""}.chk {chk_filename}{"_$SGE_TASK_ID" if multiple else ""}.fchk')

    return    
    
def extract_complexation_energy(log_file_path):
    '''
    Return the complexation energy (in kcal mol-1) from a Gaussian .log file
    looks for "(corrected)" energy, buy you could use "(raw)" if you wanted (should that be an option?)
    '''
    complexation_energy = None
    with open(log_file_path, 'r') as file:
        for line in file:
            if "complexation energy" in line:
                if "(corrected)" in line:
                    parts = line.split('=')
                    try:
                        energy_str = parts[1].split()[0]  # Get the first part after '=', which should be the energy value
                        complexation_energy = float(energy_str.strip())
                    except (IndexError, ValueError):
                        print("Error parsing complexation energy from line:", line)
    return complexation_energy

def extract_final_energy(log_file_path):
    '''
    Return the SCF energy from a Gaussian .log file (in Ha)
    '''
    final_energy = None
    with open(log_file_path, 'r') as file:
        for line in file:
            if line.startswith(" Counterpoise corrected energy ="):
                parts = line.split('=')
                try:
                    final_energy = float(parts[1].strip())
                except (IndexError, ValueError):
                    # Handle the case where parsing fails
                    print("Error parsing energy from line:", line)
                    final_energy = None
    return final_energy

def extract_translation_and_rotation_coordinates(file_path):
    '''
    Read the translation and rotation coordinates from the Gaussian log file header.
    
    ARGS:
        file_path - its the file you are reading!
    RETURNS:
        coordinates - the translation coordinates of displacement (x /y / z)
        rotations   - the rotation applied (x / y  /z)
    '''
    coordinates = []
    rotations = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip().startswith('dx='):
                match = re.search(r'dx=(-?\d+\.?\d*)/dy=(-?\d+\.?\d*)/dz=(-?\d+\.?\d*)/rx=(-?\d+\.?\d*)/ry=(-?\d+\.?\d*)/rz=(-?\d+\.?\d*)', line.strip()) # Use regex to find the dx, dy, dz, rx, ry, rz values
                if match:
                    dx, dy, dz, rx, ry, rz = map(float, match.groups())
                    coordinates.append((dx, dy, dz))
                    rotations.append((rx, ry, rz))
    return coordinates, rotations

def find_log_files(directory):
    # finds .log files and NOWT ELSE
    log_files = []
    for file in os.listdir(directory):
        if file.endswith(".log"):
            log_files.append(os.path.join(directory, file))
    return log_files
    
def clean_path(filename):
    """
    Remove files matching 'filename*.gjf' where * is an integer and 'filename.sh' in the current directory.
    """
    for log_file in glob.glob(filename.split('.')[0] + '*.gjf'):
        os.remove(log_file)

    for sh_file in glob.glob(filename.split('.')[0] + '*.sh'):
        os.remove(sh_file)    
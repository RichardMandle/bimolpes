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
	
def write_gjf(frag1, 
              frag2,
              displacement,
              nproc=4, 
              vmem=4, 
              groute='#T B3LYP cc-pVTZ EmpiricalDispersion=GD3BJ counterpoise=2',
              file_name='test'):
    '''
    Tool for writing two molecular geometries (frag1, frag2) into a gaussian .gjf input file.
    
    Args:
        frag1(2)        - fragment 1 (2) XYZ information formatted as per .gjf file with fragment=XYZ
        displacement    - job name, but with displacement coordinates encoded
        nproc           - number of CPU cores to use per job
        vmem            - ammount of ram to use (GB)function
        groute          - gaussian route information; e.g. B3LYP cc-pVTZ EmpiricalDispersion=GD3BJ counterpoise=2
        file_name       - output file name of .gjf file
    '''
    
    with open(file_name + '.gjf', 'w') as f:
        f.write('%nprocshared=' + str(nproc) + '\n')
        f.write('%mem=' + str(vmem) + 'GB\n')
        f.write(groute + '\n\n') 
        f.write(displacement + '\n\n')
        f.write('0 1\n')
        
        for i, atom in enumerate(frag1 + frag2):
            f.write(atom)
            f.write('\n')
        f.write('\n\n')
    return file_name

def make_sge_job(dir_path='noname', file ='noname', vmem=4, nproc=4, startjob=0, endjob=0):
    """
    Generate a job script for running a Gaussain job on ARC (the UoL compute clusters).
    
    Args:
        filename (str): The name of the job script file (default: 'noname').
        vmem (str): The virtual memory allocation for each job in GB(default: '8').
        nproc (int): The number of processors to use for each job (default: 8).
        startjob (int): The starting index of the task array (default: 0).
        endjob (int): The ending index of the task array (default: 0).

    Returns:
        None
    """
    
    full_path = os.path.join(dir_path, file)

    # Decide if we have a task array (i.e. multiple .gjf files) 
    if startjob + endjob == 0:
        multiple = False
        
    if startjob + endjob != 0:
        multiple = True
        
    with open(full_path + '.sh', 'w') as f:
        f.write('#$ -cwd \n')
        f.write('#$ -V\n')
        f.write('#$ -l h_rt=48:00:00\n')
        f.write('#$ -l h_vmem=' + str(vmem) + 'G\n')
        f.write('#$ -pe smp ' + str(nproc) + '\n')
        f.write('#$ -l disk=5G\n')
        if multiple:
            f.write('#$ -t ' + str(startjob) + '-' + str(endjob) + '\n')  # Create a task array

        f.write('module add gaussian\n')
        f.write('export GAUSS_SCRDIR=$TMPDIR\n')
        f.write('g16 ' + file + ('_$SGE_TASK_ID' * multiple) + '.gjf\n')
        f.write('rm *core* *.sh.* \n') #if job fails, this line cleans up messy core files
        
    return    
    
def extract_complexation_energy(log_file_path):
    '''
    Return the complexation energy (in kcal mol-1) from a Gaussian .log file
    '''
    complexation_energy = None
    with open(log_file_path, 'r') as file:
        for line in file:
            if "complexation energy" in line:
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

def extract_translation_coordinates(file_path):
    '''
    Read the translation coordinates from the Gaussian log file header (We wrote these earlier for this reason)
    '''
    coordinates = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip().startswith('dx='):
                match = re.search(r'dx=(-?\d+\.?\d*)/dy=(-?\d+\.?\d*)/dz=(-?\d+\.?\d*)', line.strip()) # Use regex to find the dx, dy, dz values
                if match:
                    dx, dy, dz = map(float, match.groups())
                    coordinates.append((dx, dy, dz))
    return coordinates

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
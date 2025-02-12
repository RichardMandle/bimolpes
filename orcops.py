# functions for orca operations (orcops)

import numpy as np
import re
import os
import glob

import geoops as geo

def write_inp(args, 
              frag1, 
              frag2,
              displacement,
              file_name='test'):
    '''
    Tool for writing two molecular geometries (frag1, frag2) into an ORCA .inp input file.
    
    Args:
        args            - arguments from bimolpes.py, incl. parameters for the calculation.
        frag1(2)        - fragment 1 (2) XYZ coordinate information.
        displacement    - job name, but with displacement coordinates encoded as dx/y/z and rotation as rx/y/z
        file_name       - output file name of .inp file
        
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

    with open(file_name + '.inp', 'w') as f:
        f.write(f"! {args.route} \n")  # ORCA functional, basis, dispersion,
        f.write(f"%pal nprocs {args.cpu} end\n")
        f.write(f"%maxcore {int(args.mem) * 1000}\n")
        f.write(f"#{displacement}\n\n")
            
        f.write("* xyz 0 1\n")  # First fragment/molecule
        for atom in frag1:
            parts = atom.split()
            formatted_parts = [parts[0]] + [format_coordinate(coord) for coord in parts[1:]]
            f.write(" ".join(formatted_parts) + "\n")
        f.write("*\n")
        
        f.write("* xyz 0 1\n")  # Second fragment
        for atom in frag2:
            parts = atom.split()
            formatted_parts = [parts[0]] + [format_coordinate(coord) for coord in parts[1:]]
            f.write(" ".join(formatted_parts) + "\n")
        f.write("*\n")   
    return file_name

def make_sge_job(args, outname, startjob=0, endjob=0):
    """
    Generate a job script for running a ORCA job on ARC (the UoL compute clusters).
    
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
        
    with open(outname + '_orca.sh', 'w') as f:
        f.write('#$ -cwd \n')
        f.write('#$ -V\n')
        f.write('#$ -l h_rt=48:00:00\n')
        f.write(f'#$ -l h_vmem={args.mem}G\n')
        f.write(f'#$ -pe smp {args.cpu}\n')
        f.write(f'#$ -l disk={args.disk}G\n')
        
        if multiple:
            f.write(f'#$ -t {startjob}-{endjob}\n')  # create a task array
            
        f.write('module add orca\n')
        f.write(f'orca {file_name}{"_$SGE_TASK_ID" if multiple else ""}.inp > {file_name}{"_$SGE_TASK_ID" if multiple else ""}.out"\n')

    return    
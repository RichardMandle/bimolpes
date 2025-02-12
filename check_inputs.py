import os
import re
import numpy as np
import processing as pro
import command_parser as cp

def extract_translation_and_rotation_from_gjf(file_path):
    """
    extracts the translation vectors (dx, dy, dz) and rotation angles (rx, ry, rz) 
    from the title section of a Gaussian .gjf file.
    """
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Identify the title line (between two empty lines)
    blank_lines = [i for i, line in enumerate(lines) if line.strip() == '']
    if len(blank_lines) < 2:
        print(f"Skipping {file_path}: Improper title format")
        return None
    
    title_line = lines[blank_lines[0] + 1].strip()
    match = re.search(r'dx=(-?\d+\.?\d*)/dy=(-?\d+\.?\d*)/dz=(-?\d+\.?\d*)/rx=(-?\d+\.?\d*)/ry=(-?\d+\.?\d*)/rz=(-?\d+\.?\d*)', title_line)
    
    if match:
        dx, dy, dz, rx, ry, rz = map(float, match.groups())
        return dx, dy, dz, rx, ry, rz
    else:
        print(f"Skipping {file_path}: No valid translation or rotation vectors found")
        return None
		
def inspect_gjf_files(args):
    """
    reads all .gjf files in the specified directory, extracts translation and rotation vectors,
    and saves them in an .npz file for visualization. The energy value is set to 42 for all entries.
    
    Energy is just saved to a made up number so we can plot it using the existing tools.
    """
    directory = args.path
    output_filename = args.out
    
    gjf_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".gjf")]
    data = []
    print(f"Reading translation/rotation vectors from {len(gjf_files)} .gjf files")
    
    for i, gjf_file in enumerate(gjf_files):
        trans_rot_vec = extract_translation_and_rotation_from_gjf(gjf_file)
        if trans_rot_vec:
            data.append((*trans_rot_vec, i))  # Set energy to index for visualising.
    
    if data:
        np.savez(output_filename, data=np.array(data))
        print(f"Saved inspection data to {output_filename}")
    else:
        print("No valid data extracted.")
#functions for analysis and/or processing

import os
import numpy as np
import zipfile
import shutil
import glob
import gauops as gps
import geoops as geo

def handle_read(args):

    # master function; reads data, saves data, prints info on minima, .
    try:
        if args.reload:
            try:
                data = reload_data_from_npz(args.filename)
            except FileNotFoundError:
                print(f"Error: The file {args.filename} does not exist.")
                return
            except Exception as e:
                print(f"An error occurred while loading the file: {e}")
                return
        else:
            try:
                data = process_files(path=args.path, method=args.method)
            except Exception as e:
                print(f"An error occurred during data processing: {e}")
                return
        
        if not args.nosave and not args.reload:
            # if we are reloading, we don't want to save the same data again
            try:
                save_data_to_npz(args.filename, data)
            except Exception as e:
                print(f"An error occurred while saving the data: {e}")
                return
            
        dx_values, dy_values, dz_values, rx_values, ry_values, rz_values, dyz_values, e_values, de_values = shape_data(data)
        
        minima_list = find_local_minima(args, data, num_minima = args.minima, e_threshold = args.ethr , d_threshold = args.dthr)
        
        minima_indices = find_indices_of_minima(data, minima_list)
        
        print_minima_report(args,minima_list, minima_indices)
        
        if args.write_minima:
            # TO DO - this code is mostly implemented, but needs testing. It needs the directory namign to be better!!!
            geo.write_minima_files(args, minima_list)

    except Exception as e:
        print(f"An unexpected error occurred: {e}")

def backup_existing_directory(directory):
    '''
    Rather than blindly overwriting existing outputs, this function will
    automatically backup a directory of grid outputs to a zip file.
    
    Args:
        directory - Directory to backup
    '''
    backup_num = 1
    backup_zip = f"{directory}_bak_{backup_num}.zip"
    
    while os.path.exists(backup_zip): # Increment the backup number until we find an available name
        backup_num += 1
        backup_zip = f"{directory}_bak_{backup_num}.zip"
    
    shutil.make_archive(directory + f"_bak_{backup_num}", 'zip', directory)
    print(f"Existing directory {directory} backed up to {backup_zip}")
        
def get_output_filename(outname=None, filename=None, filename2=None):
    """
    Generates an appropriate output filename based on the inputs.
    If 'outname' is provided, it uses it directly.
    If 'outname' is None, it tries to construct a filename using 'filename' and 'filename2'.
    If all inputs are None, it defaults to "default_filename".

    Parameters:
    - outname (str): The preferred output filename.
    - filename (str): A filename that might be used to construct a new output filename.
    - filename2 (str): Another filename that might be used to construct a new output filename.

    Returns:
    - str: The decided output filename.
    """

    if outname is None and filename is None and filename2 is None:
        outname = "default_filename"
    elif outname is None:
        if filename2 is not None:
            outname = f"{filename.split('.')[0]}_{filename2.split('.')[0]}"
        elif filename is not None:
            outname = filename.split('.')[0]
        else:
            outname = "default_filename"

    print(f"Will write data to {outname}")
    
    return outname
    
def mask_values(dx_values, dy_values, dz_values, e_values, e_max):
    '''
    Logic for masking dx/y/z where e-values are above e_max
    '''
    print(f'Masking entries where energy is <= {e_max}')
    mask = e_values <= e_max
    e_values = e_values[mask]
    dx_values = dx_values[mask]
    dy_values = dy_values[mask]
    dz_values = dz_values[mask]
    return dx_values, dy_values, dz_values, e_values
    
def zip_files(args, dir_path, zip_file):
    '''
    simple function for writing all files in a given directory to a zip file.
    matches extensions (used with .sh and .gjf/.inp to bundle things for HPC use).
    
    code looks a bit ragged but its mainly getting paths and extensions correct.
    
    Args:
        args     - parsed arguments including orca flag
        dir_path - directory to zip contents of
        zip_file - the zip file to write to
    '''

    ext = (".inp", ".sh") if args.orca else (".gjf", ".sh")

    if not zip_file.endswith('.zip'):
        zip_file = zip_file.split('.')[0] + '.zip'

    counter = np.array([0, 0])  # [count of .gjf/.inp files, count of .sh files]

    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, file)
            if os.path.isfile(file_path) and file.endswith(ext):
                zipf.write(file_path, arcname=os.path.basename(file_path))
                counter += np.array([
                    1 * file.endswith(".inp" if args.orca else ".gjf"),  # Count .gjf/.inp files
                    1 * file.endswith('.sh')  # Count .sh files
                ])

    print(f"A total of {counter[0]} {'ORCA .inp' if args.orca else 'Gaussian .gjf'} files and {counter[1]} .sh file(s) were added to {zip_file}")


def process_files(path = None, method=0):
    '''
    Wrapper function - uses the tools above to process output
    
    Args:
        path   - place to get files from
        method - if 0 = use SCF energy, if =1 get complexation energy
    
    Returns:
        data - nested list of data [dx, dy, dz, energy]
    '''
    if path is None:
        path = os.getcwd()
    files = gps.find_log_files(path)
    data = []
    print(f"Reading data from {path}")
    for file in files:
        coords, rotations = gps.extract_translation_and_rotation_coordinates(file)
        if method == 0:
            energy = gps.extract_final_energy(file)
        if method == 1:
            energy = gps.extract_complexation_energy(file)

        if coords is not None and energy is not None and rotations is not None:
            for coord, rotation in zip(coords, rotations):
                data.append([coord[0], coord[1], coord[2], rotation[0], rotation[1], rotation[2], energy])
    print(f'Processed {len(data)} configurations from {len(files)} files')
    
    return data
    
def shape_data(data):
    '''
    Just shape the data returned by process_files() for plotting
    
    Args:
        data      - the data variable to use
    
    Returns:
        dx (y,z)  - x/y/z translation coordinates
        rx (y,z)  - x/y/z rotation coordinates
        e_values  - raw energy
        de_values - delta energy
    '''
    
    dx_values = np.array([row[0] for row in data])
    dy_values = np.array([row[1] for row in data])
    dz_values = np.array([row[2] for row in data])
    rx_values = np.array([row[3] for row in data])
    ry_values = np.array([row[4] for row in data])
    rz_values = np.array([row[5] for row in data])
    dyz_values = np.sqrt((np.array([row[1] for row in data])**2) + (np.array([row[2] for row in data])**2))
    e_values = np.array([row[6] for row in data])
    de_values = np.array(e_values - np.min(e_values))
    
    return dx_values, dy_values, dz_values, rx_values, ry_values, rz_values, dyz_values, e_values, de_values

def add_npz_extension(filename):
    base, ext = os.path.splitext(filename)
    if ext != '.npz':
        filename = f"{base}.npz"
    return filename

def save_data_to_npz(filename, data):
    filename = add_npz_extension(filename)
    print(f'Saving bimolecular PES data as {filename}')
    np.savez(filename, data = data)
 
def reload_data_from_npz(filename):
    filename = add_npz_extension(filename)
    print(f'Loading bimolecular PES data from {filename}')
    reloaded_data = np.load(filename, allow_pickle = True)
    data = reloaded_data['data']
    reloaded_data.close()
    return data

def find_local_minima(args, data, num_minima=10, e_threshold=10, d_threshold=2):
    """
    Tool for finding local minima by command line.
    
    Args:
        data - data file generated by process_files function
        num_minima - the number of minima to return
        e_threshold - the energy threshold to consider; ignore points outside this
        d_threshold - distance threshold for considering minima to be unique
        args.deduplicate - If True, filters out symmetric duplicates.
    
    Returns:
        minima_list - a list of minima; their positions and energies.
    
    Geometric duplicates are checked by:
    - If they differ simply by the sign of 1 or 3 translation vectors (but not 2!)
    - If one displacement is 0, then flipping the other two can still be a duplicate.
    """
    local_data = np.array(data)
    sorted_indices = np.argsort(local_data[:, 6])  # Sort by energy
    sorted_data = local_data[sorted_indices]
    
    duplicate_count = 0
    minima_list = [sorted_data[0]]  # Always keep the global minimum

    for point in sorted_data[1:]:
        if point[3] <= (minima_list[0][3] + e_threshold) and len(minima_list) < num_minima:
            add_point = True
            for minima in minima_list:
                distance = np.sqrt(np.sum((point[:3] - minima[:3]) ** 2))
                if distance < d_threshold:
                    add_point = False
                    break
                
                if args.deduplicate:
                    # Count number of coordinate sign flips
                    sign_differences = np.sum(np.sign(point[:3]) != np.sign(minima[:3]))

                    # Check if any displacement is zero
                    zero_displacement = np.any(np.abs(point[:3]) < 1e-5)

                    # Standard duplicate check: 1 or 3 sign flips
                    is_symmetric = sign_differences in [1, 3]

                    # Special case: If one displacement is zero, allow 2 flips
                    if zero_displacement and sign_differences == 2:
                        is_symmetric = True

                    same_energy = abs(point[6] - minima[6]) < 0.01  # Energy threshold

                    if is_symmetric and same_energy:
                        add_point = False
                        duplicate_count += 1
                        break    
                        
            if add_point:
                minima_list.append(point)

        if len(minima_list) >= num_minima:
            print(f'Found maximum permitted number of local minima ({num_minima})')
            break
        
    if args.deduplicate:
        print(f'\nI found {duplicate_count} duplicate geometries when constructing the minima list.')
            
    return np.array(minima_list)

def find_indices_of_minima(data, minima):
    '''
    Given a list of data (from shape_data) and minima (from find_local_minima), return the index of these
    
    Args:
        data    -   data; taken from shape_data
        minima  -   list of local minima
    Returns:
        indices -   the indices of each minima in data.
    '''
    indices = []

    for minima_row in minima:
        match = np.where((np.array(data)[:, :3] == minima_row[:3]).all(axis=1))[0]
        if match.size > 0:
            indices.append(match[0])
        else:
            indices.append(-1)
    return indices
    
def print_minima_report(args, minima_list, minima_indices):
    '''
    Simply prints the minima_list and minima_indices passed as arguments to the terminal as a formatted table
    and writes the same to a file.
    '''
    output_filename = f'{args.filename}_minima_report.txt'
    
    with open(output_filename, 'w') as f:
        header = "\n{:<12} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7}".format(
            'Type', 'X pos.', 'Y pos.', 'Z pos.', 'rot. X', 'rot. Y', 'rot. Z', 'E_complex', 'Index')
        
        f.write(header + '\n')
        print(header)
        
        for i, minima in enumerate(minima_list):
            if i == 0: # handle the global min slightly differently
                line = "{:<12} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7}".format(
                    'Global Min.', minima[0], minima[1], minima[2], minima[3], minima[4], minima[5], minima[6], minima_indices[i])
            else:
                line = "{:<12} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7} {:<7}".format(
                    'Local Min.', minima[0], minima[1], minima[2], minima[3], minima[4], minima[5], minima[6], minima_indices[i])
            
            f.write(line + '\n')
            print(line)
    print(f'Wrote minima report to {output_filename}')
    
    
def clean_path(filename):
    """
    Remove files matching 'filename*.inp', 'filename*.gjf' where * is an integer and 'filename.sh' in the current directory.
    """
    for log_file in glob.glob(filename.split('.')[0] + '*.inp'):
        os.remove(log_file)
        
    for log_file in glob.glob(filename.split('.')[0] + '*.gjf'):
        os.remove(log_file)

    for sh_file in glob.glob(filename.split('.')[0] + '*.sh'):
        os.remove(sh_file)    
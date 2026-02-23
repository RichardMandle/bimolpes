"""
minima_analysis.py

This module provides functions to extract the final geometry and complexation energy from
Gaussian .log files, compare geometries using RMSD (with the Kabsch algorithm), and write
unique geometries (those that are not duplicates within a given RMSD threshold) to XYZ files.
It uses RDKit for obtaining atomic symbols and provides plotting functions for visualization.
"""

import os
import glob
import re
import numpy as np
from rdkit import Chem

# Get RDKit's periodic table for atomic symbols.
pt = Chem.GetPeriodicTable()


def extract_final_geometry(log_file):
    """
    Extracts the final geometry (atoms and coordinates) from a Gaussian log file.
    
    Parameters:
        log_file (str): Path to the Gaussian log file.
        
    Returns:
        atoms (list of str): List of atomic symbols.
        coords (np.ndarray): Array of shape (N, 3) with coordinates.
        
    Raises:
        ValueError: If no "Standard orientation:" block is found.
    """
    with open(log_file, 'r') as f:
        lines = f.readlines()

    final_geom_block = None
    # Loop over lines to find every "Standard orientation:" block; keep the last one.
    for i, line in enumerate(lines):
        if "Standard orientation:" in line:
            geom_block = []
            j = i + 5  # skip header lines (e.g., "Center     Atomic     Atomic...")
            while j < len(lines) and "---------------------------------------------------------------------" not in lines[j]:
                geom_block.append(lines[j].strip())
                j += 1
            final_geom_block = geom_block  # overwrite; final occurrence remains

    if final_geom_block is None:
        raise ValueError(f"No 'Standard orientation:' block found in {log_file}")

    atoms = []
    coords = []
    for line in final_geom_block:
        parts = line.split()
        if len(parts) >= 6:
            atomic_number = int(parts[1])
            # Use RDKit's periodic table to get the element symbol.
            atoms.append(pt.GetElementSymbol(atomic_number))
            coords.append([float(parts[3]), float(parts[4]), float(parts[5])])
    return atoms, np.array(coords)


def extract_final_energy(log_file):
    """
    Extracts the final complexation energy (in kcal/mol) from a Gaussian log file.
    Looks for a line containing both "complexation energy =" and "(corrected)".
    
    Parameters:
        log_file (str): Path to the Gaussian log file.
        
    Returns:
        energy (float or None): The extracted energy value, or None if not found.
    """
    energy = None
    with open(log_file, 'r') as f:
        for line in f:
            if "complexation energy =" in line and "(corrected)" in line:
                match = re.search(r'complexation energy =\s*([-+]?\d*\.\d+|\d+)', line)
                if match:
                    energy = float(match.group(1))
    return energy


def kabsch(P, Q):
    """
    Computes the optimal rotation matrix U to align set P onto set Q using the Kabsch algorithm.
    
    Parameters:
        P (np.ndarray): Array of shape (N, 3) for structure P.
        Q (np.ndarray): Array of shape (N, 3) for structure Q.
        
    Returns:
        U (np.ndarray): Rotation matrix of shape (3, 3).
    """
    # Center both sets of points.
    P_centered = P - np.mean(P, axis=0)
    Q_centered = Q - np.mean(Q, axis=0)
    
    # Covariance matrix.
    C = np.dot(P_centered.T, Q_centered)
    V, S, Wt = np.linalg.svd(C)
    
    # Correct for improper rotation (reflection) if necessary.
    d = np.linalg.det(V) * np.linalg.det(Wt)
    if d < 0:
        V[:, -1] *= -1
    U = np.dot(V, Wt)
    return U


def rmsd(P, Q):
    """
    Calculates the RMSD between two sets of points after optimal alignment.
    Assumes that the atom order in P and Q is the same.
    
    Parameters:
        P (np.ndarray): Coordinates of structure P.
        Q (np.ndarray): Coordinates of structure Q.
        
    Returns:
        float: RMSD value.
    """
    U = kabsch(P, Q)
    P_aligned = np.dot(P - np.mean(P, axis=0), U)
    Q_centered = Q - np.mean(Q, axis=0)
    diff = P_aligned - Q_centered
    return np.sqrt(np.sum(diff ** 2) / P.shape[0])


def write_xyz(filename, atoms, coords, comment):
    """
    Writes an XYZ file.
    
    Parameters:
        filename (str): Output XYZ file name.
        atoms (list of str): List of atomic symbols.
        coords (np.ndarray): Array of shape (N, 3) with coordinates.
        comment (str): Comment line (typically including a label and energy).
    """
    n_atoms = len(atoms)
    with open(filename, 'w') as f:
        f.write(f"{n_atoms}\n")
        f.write(comment + "\n")
        for atom, (x, y, z) in zip(atoms, coords):
            f.write(f"{atom} {x:.6f} {y:.6f} {z:.6f}\n")


def minima_analysis(directories, rmsd_threshold=0.2, print_mode="default", output_xyz_folder="unique_xyz"):
    """
    Processes Gaussian .log files in the given directories:
      - Extracts final geometries and complexation energies.
      - Performs pairwise RMSD comparisons (only within each directory) to detect duplicate minima.
      - Only unique geometries (those with RMSD >= threshold relative to all previous unique geometries)
        have their energy saved and are written out as XYZ files.
      
    Parameters:
        directories (list of str): List of directory paths containing Gaussian .log files.
        rmsd_threshold (float): RMSD threshold (in Angstroms) for uniqueness.
        print_mode (str): Controls printed output. Options are:
                          "verbose" - print all pairwise comparisons.
                          "default" - print only when a duplicate (RMSD < threshold) is detected.
                          "summary" - only print a summary of unique geometries.
                          "silent"  - print nothing.
        output_xyz_folder (str): Folder where unique geometry XYZ files will be saved.
        
    Returns:
        dict: A dictionary mapping each directory to a list of unique geometry records.
              Each record is a dict with keys: 'filename', 'atoms', 'coords', 'energy'.
    """
    unique_results = {}

    # Create the output folder for XYZ files if it doesn't exist.
    if not os.path.exists(output_xyz_folder):
        os.makedirs(output_xyz_folder)

    for directory in directories:
        if print_mode in ("verbose", "default"):
            print(f"\nProcessing directory: {directory}")
        log_files = glob.glob(os.path.join(directory, "*.log"))
        unique_geometries = []  # List to hold unique geometry records.
        
        for log_file in log_files:
            try:
                atoms, coords = extract_final_geometry(log_file)
                energy = extract_final_energy(log_file)
                if energy is None:
                    if print_mode in ("verbose", "default"):
                        print(f"  Warning: Energy not found in {os.path.basename(log_file)}. Skipping file.")
                    continue
            except Exception as e:
                if print_mode != "silent":
                    print(f"  Error processing {log_file}: {e}")
                continue

            duplicate_found = False
            # Compare with each unique geometry already recorded.
            for record in unique_geometries:
                current_rmsd = rmsd(coords, record["coords"])
                if print_mode == "verbose":
                    print(f"    Comparing {os.path.basename(log_file)} with {os.path.basename(record['filename'])}: RMSD = {current_rmsd:.3f} Å")
                if current_rmsd < rmsd_threshold:
                    duplicate_found = True
                    if print_mode in ("verbose", "default"):
                        print(f"    {os.path.basename(log_file)} is duplicate of {os.path.basename(record['filename'])} (RMSD = {current_rmsd:.3f} Å)")
                    break

            if not duplicate_found:
                # This geometry is unique; save its data.
                record = {
                    "filename": log_file,
                    "atoms": atoms,
                    "coords": coords,
                    "energy": energy
                }
                unique_geometries.append(record)
                if print_mode == "verbose":
                    print(f"    {os.path.basename(log_file)} is unique (Energy = {energy:.2f} kcal/mol).")
                # Write out an XYZ file for this unique geometry.
                dir_label = os.path.basename(os.path.normpath(directory))
                file_id = os.path.splitext(os.path.basename(log_file))[0]
                xyz_filename = os.path.join(output_xyz_folder, f"{dir_label}_{file_id}.xyz")
                comment = f"{dir_label} {file_id} Energy: {energy:.2f} kcal/mol"
                write_xyz(xyz_filename, atoms, coords, comment)

        unique_results[directory] = unique_geometries
        if print_mode == "summary":
            print(f"Directory: {directory} -> {len(unique_geometries)} unique geometries found.")
            for ug in unique_geometries:
                print(f"filename: {os.path.basename(ug['filename'])} with energy: {ug['energy']:.2f} kcal/mol")
                
        elif print_mode in ("verbose", "default") and len(unique_geometries) > 0:
            print(f"Summary for {directory}: {len(unique_geometries)} unique geometries.")

    return unique_results


def plot_boxplot_results(unique_results, title="Counterpoise Corrected Complexation Energy Distribution",
                           xlabel="Configuration", ylabel="Energy (kcal/mol)"):
    """
    Plots a boxplot of the complexation energies for each directory from the unique results.
    
    Parameters:
        unique_results (dict): Dictionary from minima_analysis mapping directory -> list of records.
        title (str): Plot title.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
    """
    import matplotlib.pyplot as plt

    # Build a dictionary of energies from unique_results.
    data = {}
    for directory, records in unique_results.items():
        energies = [record["energy"] for record in records if record.get("energy") is not None]
        data[directory] = energies

    # Sort keys and transform labels (e.g., take last part of path, remove "spe_" prefix)
    sorted_dirs = sorted(data.keys())
    labels = [os.path.basename(os.path.normpath(d)).split("spe_")[-1] for d in sorted_dirs]
    boxplot_data = [data[d] for d in sorted_dirs]

    plt.figure(figsize=(4, 4))
    plt.boxplot(boxplot_data, labels=labels, patch_artist=True)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.show()


def plot_swarm_results(unique_results, title="Complexation Energies for Unique Geometries",
                       xlabel="Configuration", ylabel="Energy (kcal/mol)"):
    """
    plot a swarm (strip) plot of the complexation energies for each directory from the unique results.
    
    Parameters:
        unique_results (dict): Dictionary from minima_analysis mapping directory -> list of records.
        title (str): Plot title.
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd

    # Prepare data in long form.
    data_list = []
    for directory, records in unique_results.items():
        label = os.path.basename(os.path.normpath(directory)).split("spe_")[-1]
        for record in records:
            energy = record.get("energy")
            if energy is not None:
                data_list.append((label, energy))
    df = pd.DataFrame(data_list, columns=["Configuration", "Energy"])

    plt.figure(figsize=(10, 6))
    ax = sns.swarmplot(x="Configuration", y="Energy", data=df)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()


# If run as a script, one might call minima_analysis and then plot the results.
if __name__ == "__main__":
    ## Example directories (update these paths as needed).
    #directories = [
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_rm734_anti_min",
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_rm734_par_min",
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_rm734cn_anti_min",
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_rm734cn_par_min"
    #]
    #### Run the minima analysis.
    #results = minima_analysis(directories, rmsd_threshold=0.2, print_mode="summary")
    #### Plot the results.
    #plot_boxplot_results(results)
    #plot_swarm_results(results)
    #
    ## Example directories (update these paths as needed).
    #directories = [
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_dio_anti_min",
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_dio_par_min",
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_cio_anti_min",
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_cio_par_min",
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_duzpu3f_anti_min",
    #    r"C:\RJM\bimolpes\bPES_paper\minima\min_duzpu3f_par_min"
    #]
    ## Run the minima analysis.
    #results = minima_analysis(directories, rmsd_threshold=0.05, print_mode="summary")
    ## Plot the results.
    ##plot_boxplot_results(results)
    ##plot_swarm_results(results)
    
    directories = [
        r"C:\RJM\bimolpes\bPES_paper\dio_rm734_par",
        r"C:\RJM\bimolpes\bPES_paper\dio_rm734_anti",
    ]
    # Run the minima analysis.
    results = minima_analysis(directories, rmsd_threshold=0.05, print_mode="summary")
    # Plot the results.
    plot_boxplot_results(results)
    plot_swarm_results(results)
    
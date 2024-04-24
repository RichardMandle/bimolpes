# bimolpes
Python tool for generating and analysing rigid bimolecular potential energy surfaces.

# You'll need:
Numpy, Matplotlib, scipy, Gaussian G16, access to a HPC

# What it does, and how it does it:
Bimolpes works in either `write` or `read` mode. 

With `python bimolpes.py write`:
This takes an initial geometry (molecule #1; Gaussian .log file specified as -inp) which is optimised at some level, and generates a user specified set of translation coordinates (-x -y or -z: either supply two values, min:max, or one, -val:val) The resolution is controlled with the -res flag; the analysis code is written in such a way that you can combine an initial low resolution grid with higher resolution scans over volumes of interest.

A second copy of the molecule (molecule #2; gaussian .log file is given by -inp2 or is the same as #1 by default) is imposed at each translation coordinate, provided it doesn't fall foul of our minimum and maximum atomic seperation cutoffs (min_dist, max_dist). 

You can apply rotation of molecule 2 relative to #1 using the -xa / -ya / -za flags and passing the desired angle in degrees.

To reduce the computational workload you cna use the -min and -max flags to specify the minimum and maximum atomic seperations which are permitted; those outside this range are rejected.

Each set of coordinates is then written to a Gaussian .gjf file; these include fragment information for each molecule, permitting counterpoise correction. Empirical dispersion is a essential (GD3BJ). You should use the same functional/basis set as the initial optimisation. The program writes these to the directory given by -out (or defaults to -inp1_inp2) and will write to .zip file (contains .gjf and .sh files) unless told not to.

With `python bimolpes.py read`:
After executing all .gjf files on the HPC, download these locally and store them _somewhere_. Use the -path flag to give this location. The -method flag allows to choose between SCF energy (=0) or counterpoise corrected complexation energy (=1; default).

Various options allow control of the plotted PES; -mirror mirrors the PES about the xy plane; -noplt turns off plotting; -plt_emax allows us to specify an upper limit for \Delta Energy; -plt_size controls the size of the points on the PES; -mol allows us to pass the .log file of a molecule to draw on the PES; -mol2 as before, but used for drawing a translated molecule in conjunction with -mol2_idx, which allows us to specify the index of the translation coordinates (these are outptu automatically by read mode); -real_size uses _pseudo_ realistic atom sizes based on VdW radii; -mol_size allows the base atom size to be supplied; -mol_alpha controls atom transparency; -greyscale uses a greyscale colouring of atoms; -save writes data to .npz for easy reloading; -reload loads a .npz file which is quicker than reading loads of .log files; -filename allows us to specify a filename to save/reload; -minima controls the numbeer of minima that the program will find; -ethr controls the energy cutoff for identifying discrete minima; -dthr controls the distance cutoff for identifying discrete minima.


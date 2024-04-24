# bimolpes
Python tool for generating and analysing rigid bimolecular potential energy surfaces.

# You'll need:
Numpy, Matplotlib, scipy, Gaussian G16, access to a HPC<br>

# What it does, and how it does it:
Bimolpes works in either `write` or `read` mode. <br>

With `python bimolpes.py write`:<br>
This takes an initial geometry (molecule #1; Gaussian .log file specified as `-inp`) which is optimised at some level, and generates a user specified set of translation coordinates (`-x` `-y` or `-z`: either supply two values, min:max, or one, -val:val) The resolution is controlled with the `-res` flag; the analysis code is written in such a way that you can combine an initial low resolution grid with higher resolution scans over volumes of interest.<br><br>

A second copy of the molecule (molecule #2; gaussian .log file is given by `-inp2` or is the same as #1 by default) is imposed at each translation coordinate, provided it doesn't fall foul of our minimum and maximum atomic seperation cutoffs (min_dist, max_dist). <br><br>

You can apply rotation of molecule 2 relative to #1 using the `-xa` / `-ya` / `-za` flags and passing the desired angle in degrees.<br><br>

To reduce the computational workload you cna use the `-min` and `-max` flags to specify the minimum and maximum atomic seperations which are permitted; those outside this range are rejected.<br><br>

Each set of coordinates is then written to a Gaussian .gjf file; these include fragment information for each molecule, permitting counterpoise correction. Empirical dispersion is a essential (GD3BJ). You should use the same functional/basis set as the initial optimisation. The program writes these to the directory given by `-out` (or defaults to -inp1_inp2) and will write to .zip file (contains .gjf and .sh files) unless told not to.<br><br>

With `python bimolpes.py read`:
After executing all .gjf files on the HPC, download these locally and store them _somewhere_. Use the `-path` flag to give this location. 
The `-method` flag allows to choose between SCF energy (=0) or counterpoise corrected complexation energy (=1; default).<br><br>

Various options allow control of the plotted PES; 
`-mirror` mirrors the PES about the xy plane;<br>
`-noplt` turns off plotting; <br>
`-plt_emax` allows us to specify an upper limit for \Delta Energy; <br>
`-plt_size` controls the size of the points on the PES; <br>
`-mol` allows us to pass the .log file of a molecule to draw on the PES; <br>
`-mol2` as before, but used for drawing a translated molecule in conjunction with... <br>
...`-mol2_idx`, which allows us to specify the index of the translation coordinates (these are outptu automatically by read mode); <br>
`-real_size` uses _pseudo_ realistic atom sizes based on VdW radii; <br>
`-mol_size` allows the base atom size to be supplied; <br>
`-mol_alpha` controls atom transparency; <br>
`-greyscale` uses a greyscale colouring of atoms; <br>
`-save` writes data to .npz for easy reloading; <br>
`-reload` loads a .npz file which is quicker than reading loads of .log files; <br>
`-filename` allows us to specify a filename to save/reload; <br>
`-minima` controls the numbeer of minima that the program will find... <br>
...`-ethr` controls the energy cutoff for identifying discrete minima...; <br>
...`-dthr` controls the distance cutoff for identifying discrete minima.<br>

# Worked Example:
Generate a grid of data:
`python bimolpes.py write -inp A.log -inp2 B.log -x 5 -y 5 -z 5`
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/7c09c396-cb8c-494b-b082-4a4088dc8097)
<br><br>
Execute .gjf files; use read-mode to process data and save without plotting.
`python bimolpes.py read -minima 10 -filename test -save -noplt`
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/e63623e6-1608-4f0c-8f1a-358f4a20c92b)
<br><br>
Reload the data and visualise; showing molecules 'A' and 'B' on the PES, separated according to the translation vector of the global minimum (idx = 2970)
`python bimolpes.py read -filename test -reload -mol A.log -mol2 B.log -mol2_idx 2970`
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/4abfc1bc-de31-48a2-8f54-f520a3b736e5)

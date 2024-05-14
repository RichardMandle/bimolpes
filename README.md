# bimolpes
Python tool for generating and analysing rigid bimolecular potential energy surfaces.

# You'll need:
Numpy, Mayavi, scipy, Gaussian G16, access to a HPC<br>

# What it does, and how it does it:
Bimolpes works in either `write`, `read` or `plot` mode. Argument parsing is done by the BimolPESParser class in the command_parser.py module. Some default values are stored (and can be edited) in `config.ini`<br><br>

With `python bimolpes.py write`:<br>
This takes an initial geometry (molecule #1; Gaussian .log file specified as `-inp`) which is optimised at some level, and generates a user specified set of translation coordinates (`-x` `-y` or `-z`: either supply two values, min:max, or one, -val:val) The resolution is controlled with the `-res` flag; you can supply 1 (x=y=z), 2 (x, y=z) or 3 (x, y, z) values here. The analysis code is written in such a way that you can combine an initial low resolution grid with higher resolution scans over volumes of interest, so a first pass with `-res` =1 or =2 is fine .<br><br>

If you don't really know what spacings you want then you can pass the `-est` flag; this reads the geometry in your input file (`-inp`) and finds the min/max values for x/y/z; the displacement used is 2x these x/y/z values plus the Van derWaals radii according to the atom type in question. Additionally, the user can supply an additional displacement via `-est_disp` which takes 1 (x=y=z) 2 (x, y=z) or 3 (x, y, z) values. By default, only one hemisphere is calculated, but you can calculate both by calling `-est_hemi`.<br><br>

A second copy of the molecule (molecule #2; gaussian .log file is given by `-inp2` or is the same as #1 by default) is imposed at each translation coordinate, provided it doesn't fall foul of our minimum and maximum atomic seperation cutoffs (min_dist, max_dist). <br><br>

You can apply rotation of molecule 2 relative to #1 using the `-rot` flag and passing the desired angle in degrees; again, either give 1 (x=y=z), 2 (x, y=z) or 3 (x, y, z) values here.<br><br>

To reduce the computational workload you cna use the `-min_dist` and `-max_dist` flags to specify the minimum and maximum atomic seperations which are permitted; those outside this range are rejected.<br><br>

Each set of coordinates is then written to a Gaussian .gjf file; these include fragment information for each molecule, permitting counterpoise correction. Empirical dispersion is a essential (GD3BJ). You should use the same functional/basis set as the initial optimisation. The program writes these to the directory given by `-out` (or defaults to -inp1_inp2) and will write to .zip file (contains .gjf and .sh files) unless told not to by passing the `-zip` flag. If previous grid files are found in the same output directory then these are backed up to a .zip file, unless the `-backoff` flag is called. Other options to specify here:<br>
`-cpu` = number of CPU cores (int, defaults to 4);<br>
`-mem` = amount of RAM to use (in GB; int, defaults to 4);<br>
`-groute` = Gaussian route section to use (str, defaults to #T B3LYP cc-pVTZ EmpiricalDispersion=GD3BJ counterpoise=2);<br>
`-gver` = the version of Gaussian to use (str; for example g09, g16);<br>
`-disk` = the maximum disk to use, in GB. Defaults to 5;<br>
<br><br>

With `python bimolpes.py read`:
After executing all .gjf files on the HPC, download these locally and store them _somewhere_. Use the `-path` flag to give this location. Various options are available:<br>
`-method` flag allows to choose between SCF energy (=0) or counterpoise corrected complexation energy (=1; default);<br>
`-nosave` turns off saving the read data to .npz if called; <br> 
`-filename` allows us to specify a filename to save/reload; <br>
`-reload` loads a .npz file which is quicker than reading loads of .log files; <br>
`-filename` allows us to specify a filename to save/reload (used for plotting); <br>
`-minima` controls the numbeer of minima that the program will find... <br>
...`-ethr` controls the energy cutoff for identifying discrete minima...; <br>
...`-dthr` controls the distance cutoff for identifying discrete minima.<br>
<br>
Minima information is printed to the terminal and is used to guide visualisation as needed. <br><br>
Various options allow control of the plotted PES. With `python bimolpes.py plot`:<br>
`-filename` allows us to specify a filename to reload parsed data from. <br>
<br> Various options allow control over the plotting of the PES landscape:<br>
`-plt_emax` allows us to specify an upper limit for \Delta Energy; <br>
`-plt_size` controls the size of the points on the PES (for `points` and `fancymesh` only); <br>
`-plt_alpha` controls the alpha (transparency) of points on the PES; <br>
`-plt_style` Allows choosing between normal different triangulated surface plots available through mayavi: e.g. `surface`, `points`, `wireframe`, `fancymesh` etc.; see mlab docs;<br>
`-plt_line` Set the line width used in `wireframe` or `fancymesh` plots; defaults to 2.0; set to 0 to hide lines; <br>
`-plt_cmap` name of the matplotlib cmap to use. Defaults to plasma; consult https://matplotlib.org/stable/users/explain/colors/colormaps.html;<br>
`-plt_res` controls the resolution of the resulting mlab plot object; <br>
<br> We can draw one (or two) molecules atom the PES landscape: the second molecule is displaced according to the read translation coordinates<br>
`-mol` allows us to pass the .log file of a molecule to draw on the PES; <br>
`-mol2` as before, but used for drawing a translated molecule in conjunction with... <br>
...`-mol2_idx`, which allows us to specify the index of the translation coordinates (these are outptu automatically by read mode); <br>
`-mol_real_size` uses _pseudo_ realistic atom sizes based on VdW radii; <br>
`-mol_size` allows the base atom size to be supplied; <br>
`-mol_alpha` controls atom transparency; <br>
`-mol_grey` uses a greyscale colouring of atoms; <br>

# Worked Example:
Generate a grid of data:
`python bimolpes.py write -inp A.log -inp2 B.log -x 5 -y 5 -z 5`
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/7c09c396-cb8c-494b-b082-4a4088dc8097)
<br><br>
Execute .gjf files; use read-mode to process data and save without plotting.
`python bimolpes.py read -minima 10 -filename test -save -noplt`
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/e63623e6-1608-4f0c-8f1a-358f4a20c92b)
<br><br>
Reload the data and visualise; showing molecule 'A' on the PES
`python bimolpes.py plot -filename test -mol A.log `
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/cf393ef6-566d-4882-af2e-8e81012d71f8)
<br><br>
Reload the data and visualise with some custom options; showing molecule 'rm734' on the PES in greyscale. Set the PES alpha/transparency to 0.25, use the viridis colormap; hide values with lower than -5 complexation energy; plot as points.
`python bimolpes.py plot -mol rm734.log -filename test -mol_alpha 1 -plt_alpha 0.25 -plt_cmap viridis -plt_emax -5 -greyscale -plt_style points -plt_size 10 `
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/883cb89e-6279-46ac-b9a0-cd3ce2c91903)

# bimolpes
Python tool for generating and analysing rigid bimolecular potential energy surfaces.

# You'll need:
Numpy, Mayavi, scipy, Gaussian G16 (or G09), access to a HPC<br>

# What it does, and how it does it:
Bimolpes works in either `write`, `read` or `plot` mode. Argument parsing is done by the BimolPESParser class in the command_parser.py module. Some default values are stored (and can be edited) in `config.ini`<br><br>

With `python bimolpes.py write`:<br>
This takes an initial geometry (molecule #1; Gaussian .log file specified as `-mol`) which is optimised at some level, and generates a user specified set of translation coordinates (`-x` `-y` or `-z`: either supply two values, min:max, or one, -val:val) The resolution is controlled with the `-res` flag; you can supply 1 (x=y=z), 2 (x, y=z) or 3 (x, y, z) values here. The analysis code is written in such a way that you can combine an initial low resolution grid with higher resolution scans over volumes of interest, so a first pass with `-res` =1 or =2 is fine .<br><br>

If you don't really know what spacings you want then you can pass the `-est` flag; this reads the geometry in your input file (`-inp`) and finds the min/max values for x/y/z; the displacement used is 2x these x/y/z values plus the Van derWaals radii according to the atom type in question. Additionally, the user can supply an additional displacement via `-est_disp` which takes 1 (x=y=z) 2 (x, y=z) or 3 (x, y, z) values. By default, only one hemisphere is calculated, but you can calculate both by calling `-est_hemi`.<br><br>

A second copy of the molecule (molecule #2; gaussian .log file is given by `-mol2` or is the same as #1 by default) is imposed at each translation coordinate, provided it doesn't fall foul of our minimum and maximum atomic seperation cutoffs (min_dist, max_dist). <br><br>

You can apply rotation of molecule 2 relative to #1 using the `-rot` flag and passing the desired angle in degrees; again, either give 1 (x=y=z), 2 (x, y=z) or 3 (x, y, z) values here.<br><br>

To reduce the computational workload you cna use the `-min_dist` and `-max_dist` flags to specify the minimum and maximum atomic seperations which are permitted; those outside this range are rejected. In short, this creates a 'shell' of possible translations of molecule #2 around molecule #1. <br><br>

Each set of coordinates is then written to a Gaussian .gjf file; these include fragment information for each molecule, permitting counterpoise correction. Empirical dispersion is a essential (GD3BJ). You should use the same functional/basis set as the initial optimisation. The program writes these to the directory given by `-out` (or defaults to -inp1_inp2) and will write to .zip file (contains .gjf and .sh files) unless told not to by passing the `-zip` flag. If previous grid files are found in the same output directory then these are backed up to a .zip file, unless the `-backoff` flag is called. A full list of option flags which are used in writing Gaussian files is given below:<br>
`-cpu` = number of CPU cores (int, defaults to 4);<br>
`-mem` = amount of RAM to use (in GB; int, defaults to 4);<br>
`-groute` = Gaussian route section to use (str, defaults to #T B3LYP cc-pVTZ EmpiricalDispersion=GD3BJ counterpoise=2);<br>
`-gver` = the version of Gaussian to use (str; for example g09, g16);<br>
`-disk` = the maximum disk to use, in GB. Defaults to 5;<br>
`-mol` = the .log file of moleucle #1 which will be read in;<br>
`-mol2` = the .log file of moleucle #2 which will be read in; if blank, will read a second copy of #1; <br>
`-chk` = should we write the checkpoint file? IF so, defaults to same name as the .gjf, mainly used if we need to interface to multiwfn or whatever;<br>
`-frz` = freezes atom numbers specified after the flag; useful if you want to fix position of a molecule to a vector (2 atoms) or plane (3 atoms) and perform a semi-relaxed (i.e. optimisation around constraint) calculation. <br>
<br><br>
The code is set up to auto generate .sh files for HPC submission via SGE-qsub, and will utilize task arrays for multiple jobs as apropriate. There are some options for cleaining up unwanted *.sh.* files, core* files from failed Gaussian Jobs, and for converting *.chk files to *.fchk files _via_ formchk if needed.
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
`-write_minima` writes the geometry of all minima generated above with the -minima flag to new .gjf files in a new directory. Combine this with the Gaussian options above to find local minima and reoptimise (for example).<br>
<br>
Minima information is printed to the terminal and is used to guide visualisation as needed. <br><br>
Various options allow control of the plotted PES. With `python bimolpes.py plot`:<br>
`-filename` allows us to specify a filename to reload parsed data from. <br>
<br> Various options allow control over the plotting of the PES landscape:<br>
`-plt_emax` allows us to specify an upper limit for \Delta Energy; <br>
`-plt_size` controls the size of the points on the PES (for `points` and `fancymesh` only); <br>
`-plt_alpha` controls the alpha (transparency) of points on the PES; <br>
`-plt_style` Allows choosing between normal different triangulated surface plots available through mayavi: e.g. `surface`, `points`, `wireframe`, `isosurface`, `volume`, `cut_plane` etc.; see mlab docs!<br>
`-plt_line` Set the line width used in `wireframe` or `fancymesh` plots; defaults to 2.0; set to 0 to hide lines; <br>
`-plt_cmap` name of the matplotlib cmap to use. Defaults to plasma; consult https://matplotlib.org/stable/users/explain/colors/colormaps.html;<br>
`-plt_res` controls the resolution of the resulting mlab plot object; <br>
`-plt_flipz` flip the data about the z- axis; effectively puts the molecule "on top" of the PES, a visual trick only; <br>
`-plt_minz` Only show the z-coordinate with the lowest energy, i.e. lowest energy for each combination of x/y; <br>
`-plt_axes` draw the box axes as an object; <br>
`-plt_minpoint` specify the midpoint of the colourmap; for negative values use quotation marks (e.g. "-0.05"); <br>
`-plt_min` minimum value of energy to for colourmap; for negative values use quotation marks (e.g. "-15"); <br>
`-plt_max` maximum value of energy to for colourmap; for negative values use quotation marks (e.g. "-1"); <br>
<br> We can draw one (or two) molecules atom the PES landscape: the second molecule is displaced according to the read translation coordinates<br>
`-mol` allows us to pass the .log file of a molecule to draw on the PES; <br>
`-mol2` as before, but used for drawing a translated molecule in conjunction with... <br>
...`-mol2_idx`, which allows us to specify the index of the translation coordinates (these are outptu automatically by read mode); <br>
`-mol_real_size` uses _pseudo_ realistic atom sizes based on VdW radii; <br>
`-mol_size` allows the base atom size to be supplied; <br>
`-mol_alpha` controls atom transparency; <br>
`-mol_grey` uses a greyscale colouring of atoms; <br>
`-mol_grey` colour all atoms white; <br>

<br><br>
Some more specific options for odd types of plots:
`-vis_plt` Make a vis_plot: basically, just show the extent of the translation vectors while rendering molecules #1 and #2 ontop of the grid (with the `-mol`, `-mol2` and `-mol_idx` flags <br>
`-vis_plt_x`, `-vis_plt_y`, `-vis_plt_z` control the extent in X/Y/Z; `-vis_plt_size` controls the size of the points. <br><br>

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

<br><br> Reload data, visualise the molecule RM2929 on the PES using the `-cut_plane` mode to visualise interaction strenght on a slice through the PES.
`python ..\..\bimolpes.py plot -filename ..\RM2929_par.npz -mol ..\..\CBFs\RM2931.log -plt_alpha 1 -plt_style cut_plane -plt_cmap magma`
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/e0e2eb7d-920c-45a6-a85c-556433d9f8dc)

<br><br> Generate data for the path "RM2930_anti_result", save to "RM2930_anti":
`python ..\bimolpes.py read -path RM2930_anti_result -filename RM2930_anti`
<br>
Plot the data using a midpoint of zero, min/max of -10/2 (note tue quotes around -10!), bwr cmap, draw molecule on PES, show box axes, only show minimum z- value (-plt_minz)
`python ..\bimolpes.py plot -filename RM2930_anti -plt_midpoint 0 -plt_max 2 -plt_min "-10" -plt_cmap bwr -mol ..\CBFs\RM2930.log  -plt_axes -plt_minz`
![image](https://github.com/user-attachments/assets/62b1430c-bc03-457c-9ac7-30bcf2e569db)

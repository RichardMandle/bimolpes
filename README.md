# bimolpes
repo for bimolecular potential energy surfaces

# You'll need:
Numpy, Matplotlib, scipy, Gaussian G16, access to a HPC

# What it does:
Takes an initial geometry (molecule #1) and generates a user specified set of translation coordinates (-dx:dx, -dy:dy, -dz:dz or 0:dz, resolution). A second copy of the molecule (molecule #2) is imposed at each translation coordinate, provided it doesn't fall foul of our minimum and maximum atomic seperation cutoffs (min_dist, max_dist). <br><br>
Each set of coordinates is then written to a Gaussian .gjf file; these include fragment information for each molecule, permitting counterpoise correction. Empirical dispersion is a must (GD3BJ). You should use the same functional/basis set as the initial optimisation
<br><br>
Once these calculations are done (and there could be 10k of them easily) the program extracts the translation coordinates for each geometry (which are stored in the header) and the energy (or complexation energy). These are then used to plot the bimolecular potential energy surface.

# To use:
Optimise the geometry of your molecule of interest at some level of theory (and make not of that level, it isn't automatically carried over)<br>
Read geometry and make a grid <br>`create_grid(filename = 'rm734.log', dx = 20.0, dy = 7.0, dz = 5, res = 0.5, full_dz = False, min_dist = 3, max_dist = 5)`<br>
Process the data <br>`data = process_files(path=os.getcwd(),method=0)` <br>and shape it <br>`dx_values, dy_values, dz_values, dyz_values, e_values, de_values = shape_data(data)`<br>
Get plotting in 2D <br>`make_contour_plot(dx_values,dy_values,e_values)`<br>
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/a8fdb68e-4953-4b6c-ae12-945ce861d8c1)

or in 3D <br>`make_volumetric_plot(dx_values,dy_values,dz_values,e_values)`<br>
![image](https://github.com/RichardMandle/bimolpes/assets/101199234/da63566c-f654-4a61-9524-a332f17c4121)

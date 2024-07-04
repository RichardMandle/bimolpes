# plotting functions

import numpy as np
from scipy.spatial import Delaunay
from mayavi import mlab  
from rdkit import Chem

import gauops as gps
import geoops as geo
import processing as pro

def plot_data(args):
    '''
    Master function for plotting data / molecules
    '''
    try:
        data = pro.reload_data_from_npz(args.filename)
    except FileNotFoundError:
        print(f"Error: File {args.filename} not found.")
        return
        
    dx_values, dy_values, dz_values, dyz_values, e_values, de_values = pro.shape_data(data) # shape data so its ready for plotting

    #np.savetxt("dz_values.txt", data)
        
    if args.plt_emax != 0: # if requested, mask values where energy is greater than e_min by args.plt_emax
        dx_values, dy_values, dz_values, e_values = pro.mask_values(dx_values, dy_values, dz_values, e_values, args.plt_emax)
        
    fig = mlab.figure(bgcolor=(1, 1, 1), size=(800, 600))
    
    if args.plt_flipz:
        dz_values = - dz_values # flip them so the molecule is "on top" of the surface; aesthetic feature 
    
    if not args.vis_plt:
        plot_volume(args, dx_values, dy_values, dz_values, e_values, fig = fig) # plot the volumetric surface here.
    
    if args.vis_plt:
        plot_grid(args, fig = fig) # used to visualise the grid (new function)
    
    if args.mol != None:
        plot_molecule(args, filename=args.mol, fig = fig)

    if args.mol2_idx != 0:
        if args.mol2 == None:
            filename = args.mol
        if args.mol2 != None:
            filename = args.mol2
            
        plot_molecule(args, filename=filename, displacements = data[args.mol2_idx,0:3], fig = fig)
    mlab.show()   

def plot_grid(args, fig = None):
    '''
    Little function for visuaising the grid of translations used in a bimolpes calculation, with a molecule at the centre.
    
    Probably remove this from production verseion?
    '''
    args.res = np.array([2, 2, 2])
    args.x = args.vis_plt_x
    args.y = args.vis_plt_y
    args.z = args.vis_plt_z
        
    args.min_dist = 5
    args.max_dist = 7
    
    grid = geo.gen_grid(args)

    x, y, z = zip(*grid)

    mlab.points3d(x, y, z, scale_factor=0.33, opacity = 0.25)
    

def plot_molecule(args, filename, displacements=[0, 0, 0], real_size=True, greyscale=False, fig = None):
    '''
    Plot a molecule from a Gaussian calculation using Mayavi.
    Args:
        args                - the passed arguments
        filename            - the .log file you want to open
        displacements       - set of displacements to make to the moleucle in question in x/y/z-values
        fig                 - handle of the mayavi figure to plot into
        
    args of Note:
        args.mol_size       - the base size of the moleucle atoms
        args.mol_alpha      - transparancy of the molecule(s)
        args.mol_real_size  - use "real" sized atoms, based on the VdW radii from RDKIT.
        args.mol_grey       - use (a sort of) greyscale for colouring molecule

    Returns:
        A nice molecule, rendered in the same figure as your plot.
    '''
    colors = {
        'H': (1, 1, 1),        # white
        'B': (1, 0.75, 0.8),   # pink
        'C': (0.5, 0.5, 0.5),  # gray
        'N': (0, 0, 1),        # blue
        'O': (1, 0, 0),        # red
        'F': (0, 1, 0),        # green
        'S': (1, 1, 0),        # yellow
        'Cl': (0.5, 0.5, 0),   # olive
        'Se': (0.85, 0.65, 0.13), # goldenrod
        'Br': (0.65, 0.16, 0.16), # brown
        'I': (0.5, 0, 0.5)     # purple
    }
    dx, dy, dz = displacements
    geometry = geo.generate_coords(gps.get_geometries(filename)[-1])

    if fig is None:
         fig = mlab.figure(bgcolor=(1, 1, 1), size=(800, 600))
         
    for g in geometry:
        if g[0] in colors:
            color = colors[g[0]]
        if g[0] not in colors:
            color =  (1,0,0.8) # pink! why not
        if args.mol_grey:
            color = tuple([np.mean(color)]*3)
        if args.mol_white:
            color = tuple([1,1,1])
        parts = g.split()
        if len(parts) == 4:
            sz = args.mol_size
            atom, x, y, z = parts[0], float(parts[1]) + dx, float(parts[2]) + dy, float(parts[3]) + dz
            if args.mol_real_size:
                sz = (args.mol_size ** 0.75) * (Chem.GetPeriodicTable().GetRvdw(atom))
            mlab.points3d(x, y, z, scale_factor=sz, color = color, resolution=args.plt_res, transparent=True, opacity=args.mol_alpha)

def plot_volume(args, dx_values, dy_values, dz_values, e_values, fig = None):
    '''
    Create a 3D plot only showing the minimum energy for each unique (x, y) pair over all z-values.
    
    Args:
        args:   passed arguments
        dx_values, dy_values, dz_values, e_values: Coordinate and energy arrays.
        
    args of Note:
        args.plt_size:  Size of points used in plots (scatter etc.)
        args.plt_line:  Line width to use in some plots 
        args.plt_alpha: Transparency level.
        args.plt_cmap:  Colormap name; all your favourites.
        args.plt_style: either ‘surface’ or ‘wireframe’ or ‘points’ or ‘mesh’ or ‘fancymesh’. Default: surface
        args.plt_res:   Plotting resolution (of surface, points, mesh etc).
        
    Returns:
        Makes a cool plot.
    '''  
    
    if fig is None:
        fig = mlab.figure(bgcolor=(1, 1, 1), size=(800, 600))
    
    if args.plt_max == 0:
        args.plt_max = None
        
    if args.plt_min == 0:
        args.plt_min = None

    points2D = np.vstack([dx_values, dy_values]).T
    tri = Delaunay(points2D)
    mesh = mlab.triangular_mesh(dx_values, dy_values, dz_values, tri.simplices, tube_radius = args.plt_line, scale_factor = args.plt_size, representation = args.plt_style, scalars=e_values, colormap=args.plt_cmap, opacity=args.plt_alpha, resolution = args.plt_res, vmax = args.plt_max, vmin = args.plt_min)

    if args.plt_style == 'points':
        mesh.actor.property.point_size = args.plt_size
        if args.plt_size == 1: # print a note about how the default point size probably looks bad here.
            print(f'When using plot_style points you may also want to specify a value of -plt_size (using default of 1 currently)')
        
    cb = mlab.colorbar(title='$\Delta E_{complex}$ / kcal mol$^{-1}$', orientation='horizontal')
    cb.label_text_property.color = cb.title_text_property.color = (0, 0, 0)
    cb.title_text_property.vertical_justification = 'bottom'
    cb.title_text_property.justification = 'center'
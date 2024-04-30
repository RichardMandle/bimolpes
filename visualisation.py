# plotting functions

import numpy as np
from mayavi import mlab
from scipy.spatial import Delaunay

from rdkit import Chem

import gauops as gps
import geoops as geo

def plot_molecule_mayavi(filename, displacements=[0, 0, 0], size=2, alpha=1, real_size=True, greyscale=False, fig = None):
    '''
    Plot a molecule from a Gaussian calculation using Mayavi.
    Args:
        filename - the .log file you want to open
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
        if greyscale:
            color = tuple([np.mean(color)]*3)
        parts = g.split()
        if len(parts) == 4:
            sz = size
            atom, x, y, z = parts[0], float(parts[1]) + dx, float(parts[2]) + dy, float(parts[3]) + dz
            if real_size == True:
                sz = (size ** 0.75) * (Chem.GetPeriodicTable().GetRvdw(atom))
            mlab.points3d(x, y, z, scale_factor=sz, color = color, resolution=20, transparent=True, opacity=alpha)

def plot_volumetric_mayavi(dx_values, dy_values, dz_values, e_values, e_max=0, size=1, line_width = 2.0, alpha=0.5, plot_style='surface', cmap='plasma', fig = None):
    '''
    Create a 3D plot only showing the minimum energy for each unique (x, y) pair over all z-values.
    
    Args:
        dx_values, dy_values, dz_values, e_values: Coordinate and energy arrays.
        e_max (float): Maximum energy value to display.
        size (float): Size of points for scatter plot.
        mirror (bool): Whether to mirror the plot across the z=0 plane.
        alpha (float): Transparency level.
        plot_style (str): either ‘surface’ or ‘wireframe’ or ‘points’ or ‘mesh’ or ‘fancymesh’. Default: surface
        cmap (str): Colormap name.
    '''
    if fig is None:
        fig = mlab.figure(bgcolor=(1, 1, 1), size=(800, 600))

    points2D = np.vstack([dx_values, dy_values]).T
    tri = Delaunay(points2D)
    mesh = mlab.triangular_mesh(dx_values, dy_values, dz_values, tri.simplices, tube_radius = line_width, scale_factor = size, representation = plot_style, scalars=e_values, colormap=cmap, opacity=alpha, resolution = 20)

    if plot_style == 'points':
        mesh.actor.property.point_size = size
        if size == 1: # print a note about how the default point size probably looks bad here.
            print(f'When using plot_style points you may also want to specify a value of plt_size (using default of 1 currently)')
        
    cb = mlab.colorbar(title='$\Delta E_{complex}$ / kcal mol$^{-1}$', orientation='horizontal')
    cb.label_text_property.color = cb.title_text_property.color = (0, 0, 0)
    cb.title_text_property.vertical_justification = 'bottom'
    cb.title_text_property.justification = 'center'
    
def finalize_scene():
    ''' Just finalize the current mlab scene and nowt else'''
    mlab.show()
# plotting functions

import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable, plasma
from matplotlib.tri import Triangulation

from rdkit import Chem

import gauops as gps
import geoops as geo

def atomic_number_to_grayscale(z, max_z = 12):
    grayscale_value = min((z + (max_z / 2)) / max_z, 1)
    return (grayscale_value, grayscale_value, grayscale_value)

def plot_molecule(filename, displacements = [0,0,0], ax = None, size = 256, alpha = 1, real_size = True, greyscale = False):
    '''
    Roughly plot a molecule from a Gaussian calculation using matplotlib.
    Args:
        filename - the .log file you want to open
        ax       - the axis you want to draw the molecule onto.
        
    example usage:  
        plot_molecule(filename='rm734.log')
    '''
    
    colors = {'H':'white','B':'pink','C':'gray','N':'blue','O':'red','F':'green','S':'yellow','Cl':'olive', 'Se':'goldenrod','Br':'brown','I':'purple'} # add new colours here as needed...
    dx, dy, dz = displacements[0], displacements[1], displacements[2], 
    geometry = geo.generate_coords(gps.get_geometries(filename)[-1])
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
    for g in geometry:
        if g[0] in colors:
            color = colors[g[0]]
        if g[0] not in colors:
            color = 'pink' # pink! why not
        
        if greyscale == True:
            color = atomic_number_to_grayscale(Chem.GetPeriodicTable().GetAtomicNumber(g[0]))
            
        parts = g.split()
        if len(parts) == 4:
            sz = size
            atom, x, y, z = parts[0], float(parts[1]) + dx, float(parts[2]) + dy, float(parts[3]) + dz
            if real_size == True:
                sz = (size ** 0.75) * (Chem.GetPeriodicTable().GetRvdw(atom))
            ax.scatter(x, y, z, color = color, s = sz, alpha = alpha, edgecolor = 'black')

      
def plot_volumetric(dx_values, dy_values, dz_values, e_values, e_max = 0, size = 1, mirror = False, alpha = 0.5, plot_style = 'trisurf', cmap = 'plasma', plot_mol = False):
    '''
    Create a 3D plot only showing the minimum energy for each unique (x, y) pair over all z-values.
    
    example:
    ax = plot_volumetric(dx_values, dy_values, dz_values, e_values, e_max = 0, size = 1, mirror = True, alpha = 1)
    
    can plot a molecule on top using the returned axis handle:
    plot_molecule(filename='rm734.log', ax = ax, size = 128, alpha = 1)
    
    NOTE - this function is huge and doing loads of things, could be broken down a bit.
    '''
    if plot_style not in ['scatter', 'trisurf']:
        raise ValueError(f"Unsupported plot type '{plot_style}'. Supported types are 'scatter' and 'trisurf'.")

    
    data = np.core.records.fromarrays([dx_values, dy_values, dz_values, e_values],
                                      names='dx, dy, dz, e')
    
    if e_max != 0:
        data['e'][data['e'] > e_max] = np.nan
    
    unique_xy = np.unique(data[['dx', 'dy']], axis=0)
    min_energy_points = np.array([(x, y, np.min(data['dz'][(data['dx'] == x) & (data['dy'] == y)]),
                                  np.min(data['e'][(data['dx'] == x) & (data['dy'] == y)]))
                                 for x, y in unique_xy])

    dx_min, dy_min, dz_min, e_min = min_energy_points[:, 0], min_energy_points[:, 1], min_energy_points[:, 2], min_energy_points[:, 3]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    norm = Normalize(vmin=np.nanmin(e_values), vmax=np.nanmax(e_max))
    mappable = ScalarMappable(norm=norm, cmap = cmap)
    colors = mappable.to_rgba(e_values)

    if plot_style == 'scatter':
        sc = ax.scatter(dx_min, dy_min, dz_min, s=size, c=e_min.astype(float), cmap = cmap, alpha=alpha)
        if mirror == True:
            sc = ax.scatter(dx_min, dy_min, dz_min, s=size, c=e_min.astype(float), cmap = cmap, alpha=alpha)

    elif plot_style == 'trisurf':
        if (e_max != 0):
            colors[np.array(e_values) > e_max, 3] = 0 # set alpha to 0 where energy exceeds e_max; basically a hacky way of handling the trisurf grid
            colors[np.array(e_values) <= e_max, 3] = colors[np.array(e_values) <= e_max, 3] * alpha # Now pass the user supplied alpha           
        tri = Triangulation(dx_values, dy_values)
        trisurf = ax.plot_trisurf(tri, list(-np.array(dz_values)), color='white')  # Initial dummy color; we'll bin that off shortly
        trisurf.set_facecolors(colors[tri.triangles].mean(axis=1))

        if mirror:
            trisurf_mirror = ax.plot_trisurf(tri, dz_values, alpha=alpha, color='white')
            trisurf_mirror.set_facecolors(colors[tri.triangles].mean(axis=1))

    plt.colorbar(mappable, label='$\Delta E_{complex}$ / kcal mol$^{-1}$', shrink=0.45)

    x_range = np.ptp(dx_min)  # Peak to peak (range) for dx
    y_range = np.ptp(dy_min)  # Peak to peak (range) for dy
    z_range = np.ptp(dz_min)  # Peak to peak (range) for dz
    max_range = np.array([x_range, y_range, z_range]).max()
    aspect_ratio = max_range / np.array([x_range, y_range, z_range])

    ax.set_box_aspect(1/aspect_ratio)  # Set the aspect ratio based on the ranges
    ax.axis('off')
    ax.view_init(elev=90, azim=0) # this is a nice view.
    
    if plot_mol == False: 
        plt.show() # if we aren't plottin' a molecule do it here.

    return ax
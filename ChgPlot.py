from jasp import *
from enthought.mayavi import mlab
from ase.data.molecules import molecule
from ase.data import vdw_radii
from ase.data.colors import cpk_colors
import numpy as np

def ChargeDensityFog3D(directory, save_file=None , DrawAtoms=True, DrawCell=True, opacity_range=[0,1]):
    # Get data from calculation files
    with jasp(directory) as calc:
        atoms = calc.get_atoms()
        x, y, z, cd = calc.get_charge_density()

    mlab.figure(bgcolor=(0,0,0), size=(640,480))

    # Draw atoms
    if DrawAtoms == True:
        for atom in atoms:
            mlab.points3d(atom.x,
                          atom.y,
                          atom.z,
                          scale_factor=vdw_radii[atom.number]/10.,
                          resolution=16,
                          color=tuple(cpk_colors[atom.number]),
                          scale_mode='none')
    # Draw unit cell
    if DrawCell == True:
        a1, a2, a3 = atoms.get_cell()
        origin = [0,0,0]
        cell_matrix = [[origin, a1],
                       [origin, a2],
                       [origin, a3],
                       [a1, a1+a2],
                       [a1, a1+a3],
                       [a2, a2+a1],
                       [a2, a2+a3],
                       [a3, a1+a3],
                       [a3, a2+a3],
                       [a1+a2, a1+a2+a3],
                       [a2+a3, a1+a2+a3],
                       [a1+a3, a1+a3+a2]] # contains all points on the box
        for p1, p2 in cell_matrix:
            mlab.plot3d([p1[0], p2[0]], # x-coords of box
                        [p1[1], p2[1]], # y-coords
                        [p1[2], p2[2]]) # z-coords

    # Plot the charge density
    src = mlab.pipeline.scalar_field(x, y, z, cd) #Source data
    vmin = cd.min() #find minimum and maximum value of CD data
    vmax = cd.max()
    vol = mlab.pipeline.volume(src) # Make a volumetric representation of the data

    # Set opacity transfer function
    from tvtk.util.ctf import PiecewiseFunction
    otf = PiecewiseFunction()
    otf.add_point(vmin, opacity_range[0]) #Transparency at zero electron density
    otf.add_point(vmax*1, opacity_range[1]) #Transparency at max electron density
    vol._otf=otf
    vol._volume_property.set_scalar_opacity(otf)

    #Show a legend
    mlab.colorbar(title="e- density\n(e/Ang^3)", orientation="vertical", nb_labels=5, label_fmt='%.2f')
    mlab.view(azimuth=-90, elevation=90, distance='auto') # Set viewing angle
    mlab.show()
    if save_file != None:
        mlab.savefig(save_file)

def ChargeDensitySlice3D(directory, save_file=None, x_loc=None, y_loc=None, z_loc=None, DrawAtoms=True, DrawCell=True):
    # Get data from calculation files
    with jasp(directory) as calc:
        atoms = calc.get_atoms()
        x, y, z, cd = calc.get_charge_density()

    mlab.figure(bgcolor=(0,0,0), size=(640,480))

    # Draw atoms
    if DrawAtoms == True:
        for atom in atoms:
            mlab.points3d(atom.x,
                          atom.y,
                          atom.z,
                          scale_factor=vdw_radii[atom.number]/10.,
                          resolution=16,
                          color=tuple(cpk_colors[atom.number]),
                          scale_mode='none')
    # Draw unit cell
    if DrawCell == True:
        a1, a2, a3 = atoms.get_cell()
        origin = [0,0,0]
        cell_matrix = [[origin, a1],
                       [origin, a2],
                       [origin, a3],
                       [a1, a1+a2],
                       [a1, a1+a3],
                       [a2, a2+a1],
                       [a2, a2+a3],
                       [a3, a1+a3],
                       [a3, a2+a3],
                       [a1+a2, a1+a2+a3],
                       [a2+a3, a1+a2+a3],
                       [a1+a3, a1+a3+a2]] # contains all points on the box
        for p1, p2 in cell_matrix:
            mlab.plot3d([p1[0], p2[0]], # x-coords of box
                        [p1[1], p2[1]], # y-coords
                        [p1[2], p2[2]]) # z-coords

    # Plot the charge density on three perpendicular planes at the center of the cell
    src = mlab.pipeline.scalar_field(x, y, z, cd) #Source data
    vmin = cd.min() #find minimum and maximum value of CD data
    vmax = cd.max()
    if x_loc != None:
        sliceX = mlab.pipeline.image_plane_widget(src,
                                                  plane_orientation = 'x_axes',
                                                  slice_index = x_loc * 10,
                                                  transparent = True)
    if y_loc != None:
        sliceY = mlab.pipeline.image_plane_widget(src,
                                                  plane_orientation = 'y_axes',
                                                  slice_index = y_loc * 10,
                                                  transparent = True)
    if z_loc != None:
        sliceZ = mlab.pipeline.image_plane_widget(src,
                                                  plane_orientation = 'z_axes',
                                                  slice_index = z_loc * 10,
                                                  transparent = True)
    #Show a legend
    mlab.colorbar(title="e- density\n(e/Ang^3)", orientation="vertical", nb_labels=5, label_fmt='%.2f')
    mlab.show()
    if save_file != None:
        mlab.savefig(save_file)

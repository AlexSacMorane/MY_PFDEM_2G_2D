# -*- encoding=utf-8 -*-

from yade import pack, utils, plot, export
from potential_utils import *
from pathlib import Path
import math, os, errno, shutil
import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Functions

def run_moose():
    '''
    Prepare and run moose simulation.
    '''
    # from dem to pf
    move_phasefield() # in dem_to_pf.py
    compute_kc() # in dem_to_pf.py
    compute_ed() # in dem_to_pf.py
    write_i() # in dem_to_pf.py

    # pf
    os.system('mpiexec -n '+str(n_proc)+' ~/projects/moose/modules/phase_field/phase_field-opt -i pf.i')

    # from pf to dem
    last_j_str = sort_files() # in pf_to_dem.py

    print('Reading data')
    read_vtk(last_j_str) # in pf_to_dem.py

# ------------------------------------------------------------------------------------------------------------------------------------------ #

# Own functions
execfile('dem_to_pf.py')
execfile('pf_to_dem.py')
execfile('ic.py')
execfile('dem.py')
execfile('tools.py')

# Parameters
execfile('Parameters.py')

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Plan

# folders
create_folder('vtk') # from tools.py
create_folder('plot') # from tools.py
create_folder('plot/contact') # from tools.py
create_folder('data') # from tools.py
create_folder('input') # from tools.py
# vtk exporter
vtkExporter = export.VTKExporter('vtk/2grains')
vtkExporter_1 = export.VTKExporter('vtk/grain_1')
vtkExporter_2 = export.VTKExporter('vtk/grain_2')
# Plot
create_plots()
# Engines
create_engines()
# the list of the nodes
L_XYZ = None
L_L_i_XYZ_not_used = None

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Create initial condition

create_materials() # from ic.py
create_2_spheres() # from ic.py
create_solute() # from ic.py
compute_dt() # from ic.py

# Plot
fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(16,9))
# eta 1
im = ax1.imshow(eta_1_map, interpolation = 'nearest', extent=(x_L[0],x_L[-1],y_L[0],y_L[-1]))
fig.colorbar(im, ax=ax1)
ax1.set_title(r'Map of $\eta_1$',fontsize = 30)
# eta 2
im = ax2.imshow(eta_2_map, interpolation = 'nearest', extent=(x_L[0],x_L[-1],y_L[0],y_L[-1]))
fig.colorbar(im, ax=ax2)
ax2.set_title(r'Map of $\eta_2$',fontsize = 30)
# solute
im = ax3.imshow(c_map, interpolation = 'nearest', extent=(x_L[0],x_L[-1],y_L[0],y_L[-1]))
fig.colorbar(im, ax=ax3)
ax3.set_title(r'Map of solute',fontsize = 30)
# close
fig.savefig('plot/Map_etas_solute.png')
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# PFDEM iteration

i_DEMPF_ite = 0
while i_DEMPF_ite < n_DEMPF_ite:
    # DEM
    print('Running DEM')
    O.run()
    O.wait()

    # DEM->PF, PF, PF->DEM
    run_moose()

    # prepare next iteration
    pos_1 = [O.bodies[0].state.pos[0], O.bodies[0].state.pos[1], O.bodies[0].state.pos[2]]
    pos_2 = [O.bodies[1].state.pos[0], O.bodies[1].state.pos[1], O.bodies[1].state.pos[2]]

    O.reset()
    plot.resetData()
    create_materials() # from ic.py
    create_potential_block() # from dem.py
    create_engines() # from ic.py
    compute_dt() # from ic.py
    i_DEMPF_ite = i_DEMPF_ite + 1

    L_sum_eta_1.append(np.sum(eta_1_map))
    L_sum_eta_2.append(np.sum(eta_2_map))
    L_sum_c.append(np.sum(c_map))

    # plot eta_i, c
    fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
    ax1.plot(L_sum_eta_1)
    ax2.plot(L_sum_eta_2)
    ax3.plot(L_sum_c)
    fig.savefig('plot/etai_c.png')
    plt.close(fig)

    # plot displacement
    L_disp_init = [0]
    L_disp = [0]
    for i_disp in range(len(L_displacement)):
        L_disp_init.append(L_disp_init[-1]+L_displacement[i_disp])
        if i_disp >= 1:
            L_disp.append(L_disp[-1]+L_displacement[i_disp])

    fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(16,9))
    ax1.plot(L_disp)
    ax2.plot(L_disp_init)
    fig.savefig('plot/displacement.png')
    plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# close simulation

print('Simulation ends')

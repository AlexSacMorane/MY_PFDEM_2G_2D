# -*- encoding=utf-8 -*-

import math, os, errno, pickle
import numpy as np
import matplotlib.pyplot as plt

# own
from dem_to_pf import *
from pf_to_dem import *
from ic import *
from dem import *
from tools import *
from Parameters import *

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Functions

def run_moose(dict_user, dict_sample):
    '''
    Prepare and run moose simulation.
    '''
    # from dem to pf
    move_phasefield(dict_user, dict_sample) # in dem_to_pf.py
    compute_kc(dict_user, dict_sample) # in dem_to_pf.py
    compute_ed(dict_user, dict_sample) # in dem_to_pf.py
    write_i(dict_user, dict_sample) # in dem_to_pf.py

    # pf
    print('Running PF')
    os.system('mpiexec -n '+str(dict_user['n_proc'])+' ~/projects/moose/modules/phase_field/phase_field-opt -i pf.i')

    # from pf to dem
    last_j_str = sort_files(dict_user, dict_sample) # in pf_to_dem.py

    print('Reading data')
    read_vtk(dict_user, dict_sample, last_j_str) # in pf_to_dem.py

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def run_yade(dict_user, dict_sample):
    '''
    Prepare and run yade simulation.
    '''
    # from pf to dem
    compute_plane(dict_user, dict_sample) # from pf_to_dem.py
    # transmit data
    dict_save = {
    'Kn': dict_user['Kn'],
    'Ks': dict_user['Ks'],
    'force_applied': dict_user['force_applied'],
    'pos_1': dict_sample['pos_1'],
    'pos_2': dict_sample['pos_2'],
    'r_pot': dict_user['r_pot'],
    'n_ite_max': dict_user['n_ite_max'],
    'steady_state_detection': dict_user['steady_state_detection'],
    'n_steady_state_detection': dict_user['n_steady_state_detection'],
    'i_DEMPF_ite': dict_sample['i_DEMPF_ite']
    }
    with open('data/main_to_dem.data', 'wb') as handle:
        pickle.dump(dict_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # dem
    print('Running DEM')
    os.system('yadedaily -j '+str(dict_user['n_proc'])+' -x -n dem_base.py')

    # sort files
    sort_files_yade() # in dem.py

    # load data
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    dict_sample['pos_1'] = dict_save['pos_1']
    dict_sample['pos_2'] = dict_save['pos_2']

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Plan

# folders
create_folder('vtk') # from tools.py
create_folder('plot') # from tools.py
create_folder('plot/contact') # from tools.py
create_folder('data') # from tools.py
create_folder('input') # from tools.py

# get parameters
dict_user = get_parameters() # from Parameters.py
dict_sample = {}

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Create initial condition

create_2_spheres(dict_user, dict_sample) # from ic.py
create_solute(dict_user, dict_sample) # from ic.py

# compute tracker
dict_user['L_sum_eta_1'].append(np.sum(dict_sample['eta_1_map']))
dict_user['L_sum_eta_2'].append(np.sum(dict_sample['eta_2_map']))
dict_user['L_sum_c'].append(np.sum(dict_sample['c_map']))

# Plot
fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(16,9))
# eta 1
im = ax1.imshow(dict_sample['eta_1_map'], interpolation = 'nearest', extent=(dict_sample['x_L'][0],dict_sample['x_L'][-1],dict_sample['y_L'][0],dict_sample['y_L'][-1]))
fig.colorbar(im, ax=ax1)
ax1.set_title(r'Map of $\eta_1$',fontsize = 30)
# eta 2
im = ax2.imshow(dict_sample['eta_2_map'], interpolation = 'nearest', extent=(dict_sample['x_L'][0],dict_sample['x_L'][-1],dict_sample['y_L'][0],dict_sample['y_L'][-1]))
fig.colorbar(im, ax=ax2)
ax2.set_title(r'Map of $\eta_2$',fontsize = 30)
# solute
im = ax3.imshow(dict_sample['c_map'], interpolation = 'nearest', extent=(dict_sample['x_L'][0],dict_sample['x_L'][-1],dict_sample['y_L'][0],dict_sample['y_L'][-1]))
fig.colorbar(im, ax=ax3)
ax3.set_title(r'Map of solute',fontsize = 30)
# close
fig.savefig('plot/Map_etas_solute.png')
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# PFDEM iteration

dict_sample['i_DEMPF_ite'] = 0
while dict_sample['i_DEMPF_ite'] < dict_user['n_DEMPF_ite']:
    dict_sample['i_DEMPF_ite'] = dict_sample['i_DEMPF_ite'] + 1

    # DEM
    run_yade(dict_user, dict_sample)

    # DEM->PF, PF, PF->DEM
    run_moose(dict_user, dict_sample)

    # compute tracker
    dict_user['L_sum_eta_1'].append(np.sum(dict_sample['eta_1_map']))
    dict_user['L_sum_eta_2'].append(np.sum(dict_sample['eta_2_map']))
    dict_user['L_sum_c'].append(np.sum(dict_sample['c_map']))
    dict_user['L_sum_mass'].append(np.sum(dict_sample['eta_1_map'])+np.sum(dict_sample['eta_2_map'])+np.sum(dict_sample['c_map']))

    # plot eta_i, c
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))
    ax1.plot(dict_user['L_sum_eta_1'])
    ax1.set_title(r'$\Sigma\eta_1$')
    ax2.plot(dict_user['L_sum_eta_2'])
    ax2.set_title(r'$\Sigma\eta_2$')
    ax3.plot(dict_user['L_sum_c'])
    ax3.set_title(r'$\Sigma C$')
    ax4.plot(dict_user['L_sum_mass'])
    ax4.set_title(r'$\Sigma\eta_1 + \Sigma\eta_2 + \Sigma c$')
    fig.savefig('plot/etai_c.png')
    plt.close(fig)

    # plot displacement
    L_disp_init = [0]
    L_disp = [0]
    for i_disp in range(len(dict_user['L_displacement'])):
        L_disp_init.append(L_disp_init[-1]+dict_user['L_displacement'][i_disp])
        if i_disp >= 1:
            L_disp.append(L_disp[-1]+dict_user['L_displacement'][i_disp])
    fig, (ax1,ax2) = plt.subplots(nrows=1,ncols=2,figsize=(16,9))
    ax1.plot(L_disp)
    ax2.plot(L_disp_init)
    fig.savefig('plot/displacement.png')
    plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# close simulation

print('Simulation ends')

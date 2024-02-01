# -*- encoding=utf-8 -*-

import math, os, errno, pickle
import numpy as np
import matplotlib.pyplot as plt
import time

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
    tic_dem_to_pf = time.perf_counter() # compute dem_to_pf performances
    move_phasefield(dict_user, dict_sample) # in dem_to_pf.py
    compute_kc(dict_user, dict_sample) # in dem_to_pf.py
    compute_ed(dict_user, dict_sample) # in dem_to_pf.py
    write_i(dict_user, dict_sample) # in dem_to_pf.py
    tac_dem_to_pf = time.perf_counter() # compute dem_to_pf performances
    dict_user['L_t_dem_to_pf'].append(tac_dem_to_pf-tic_dem_to_pf)

    # compare contact volume in Moose and in Yade
    compare_volumes(dict_user, dict_sample) # in dem_to_pf.py

    # pf
    print('Running PF')
    tic_pf = time.perf_counter() # compute pf performances
    os.system('mpiexec -n '+str(dict_user['n_proc'])+' ~/projects/moose/modules/phase_field/phase_field-opt -i pf.i')
    tac_pf = time.perf_counter() # compute pf performances
    dict_user['L_t_pf'].append(tac_pf-tic_pf)

    # from pf to dem
    last_j_str = sort_files(dict_user, dict_sample) # in pf_to_dem.py

    print('Reading data')
    tic_pf_to_dem = time.perf_counter() # compute pf_to_dem performances
    read_vtk(dict_user, dict_sample, last_j_str) # in pf_to_dem.py
    tac_pf_to_dem = time.perf_counter() # compute pf_to_dem performances
    dict_user['L_t_pf_to_dem_2'].append(tac_pf_to_dem-tic_pf_to_dem)

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def run_yade(dict_user, dict_sample):
    '''
    Prepare and run yade simulation.
    '''
    # from pf to dem
    tic_pf_to_dem = time.perf_counter() # compute pf_to_dem performances
    compute_vertices(dict_user, dict_sample) # from pf_to_dem.py
    tac_pf_to_dem = time.perf_counter() # compute pf_to_dem performances
    dict_user['L_t_pf_to_dem_1'].append(tac_pf_to_dem-tic_pf_to_dem)

    # shape evolution
    with open('data/planes.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    # save initial shapes
    if dict_user['L_vertices_1_init'] == None:
        dict_user['L_vertices_1_init'] = dict_save['L_vertices_1']
        dict_user['L_vertices_2_init'] = dict_save['L_vertices_2']
    #compare current shape and initial one
    else :
        fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))
        # g1
        L_x, L_y = tuplet_to_list(dict_user['L_vertices_1_init']) # from tools.py
        ax1.plot(L_x, L_y, label='Initial')
        L_x, L_y = tuplet_to_list(dict_save['L_vertices_1']) # from tools.py
        ax1.plot(L_x, L_y, label='Current')
        ax1.legend()
        ax1.axis('equal')
        ax1.set_title(r'G1',fontsize=20)
        # g2
        L_x, L_y = tuplet_to_list(dict_user['L_vertices_2_init']) # from tools.py
        ax2.plot(L_x, L_y, label='Initial')
        L_x, L_y = tuplet_to_list(dict_save['L_vertices_2']) # from tools.py
        ax2.plot(L_x, L_y, label='Current')
        ax2.legend()
        ax2.axis('equal')
        ax2.set_title(r'G2',fontsize=20)
        # close
        plt.suptitle('Shapes evolution', fontsize=20)
        fig.savefig('plot/shape_evolution.png')
        plt.close(fig)

    # transmit data
    dict_save = {
    'E': dict_user['E'],
    'Poisson': dict_user['Poisson'],
    'force_applied': dict_user['force_applied'],
    'pos_1': dict_sample['pos_1'],
    'pos_2': dict_sample['pos_2'],
    'n_ite_max': dict_user['n_ite_max'],
    'steady_state_detection': dict_user['steady_state_detection'],
    'n_steady_state_detection': dict_user['n_steady_state_detection'],
    'i_DEMPF_ite': dict_sample['i_DEMPF_ite']
    }
    with open('data/main_to_dem.data', 'wb') as handle:
        pickle.dump(dict_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # dem
    print('Running DEM')
    tic_dem = time.perf_counter() # compute dem performances
    os.system('yadedaily -j '+str(dict_user['n_proc'])+' -x -n dem_base.py')
    tac_dem = time.perf_counter() # compute dem performances
    dict_user['L_t_dem'].append(tac_dem-tic_dem)

    # sort files
    sort_files_yade() # in dem_to_pf.py

    # load data
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    dict_sample['pos_1'] = dict_save['pos_1']
    dict_sample['pos_2'] = dict_save['pos_2']

    # tracker
    dict_user['L_distance_extrema'].append(dict_save['distance_extrema'])
    dict_user['L_equivalent_area'].append(dict_save['equivalent_area'])
    dict_user['L_contact_overlap'].append(dict_save['contact_overlap'])
    dict_user['L_contact_area'].append(dict_save['contact_area'])
    dict_user['L_contact_volume_yade'].append(dict_save['contact_volume'])

    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(16,9))
    # overlap
    ax1.plot(dict_user['L_distance_extrema'], label='Distance vertices')
    ax1.set_title(r'Distance between extrema ($m$)',fontsize=20)
    # surface
    ax2.plot(dict_user['L_equivalent_area'])
    ax2.set_title(r'Equivalent surface ($m^2$)',fontsize=20)
    # volume
    ax3.plot(dict_user['L_contact_volume_yade'])
    ax3.set_title(r'Volume ($m^3$)',fontsize=20)
    # close
    plt.suptitle('Contact', fontsize=20)
    fig.savefig('plot/contact_v_s_d.png')
    plt.close(fig)

    # tracker
    dict_user['L_n_v_1'].append(dict_save['n_v_1'])
    dict_user['L_n_v_2'].append(dict_save['n_v_2'])
    dict_user['L_n_v_1_target'].append(dict_save['n_v_1_target'])
    dict_user['L_n_v_2_target'].append(dict_save['n_v_2_target'])

    # plot
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    ax1.plot(dict_user['L_n_v_1'], label='N vertices g1', color='r')
    ax1.plot(dict_user['L_n_v_1_target'], label='N vertices g1 targetted', color='r', linestyle='dotted')
    ax1.plot(dict_user['L_n_v_2'], label='N vertices g2', color='b')
    ax1.plot(dict_user['L_n_v_2_target'], label='N vertices g2 targetted', color='b', linestyle='dotted')
    ax1.legend()
    fig.savefig('plot/n_vertices.png')
    plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Plan

# folders
create_folder('vtk') # from tools.py
create_folder('plot') # from tools.py
create_folder('plot/contact') # from tools.py
create_folder('data') # from tools.py
create_folder('input') # from tools.py
create_folder('dict') # from tools.py

# get parameters
dict_user = get_parameters() # from Parameters.py
dict_sample = {}

# compute performances
tic = time.perf_counter()

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Create initial condition

create_2_spheres(dict_user, dict_sample) # from ic.py
create_solute(dict_user, dict_sample) # from ic.py

# compute tracker
dict_user['L_sum_eta_1'].append(np.sum(dict_sample['eta_1_map']))
dict_user['L_sum_eta_2'].append(np.sum(dict_sample['eta_2_map']))
dict_user['L_sum_c'].append(np.sum(dict_sample['c_map']))
dict_user['L_sum_mass'].append(np.sum(dict_sample['eta_1_map'])+np.sum(dict_sample['eta_2_map'])+np.sum(dict_sample['c_map']))
dict_user['L_m_eta_1'].append(np.mean(dict_sample['eta_1_map']))
dict_user['L_m_eta_2'].append(np.mean(dict_sample['eta_2_map']))
dict_user['L_m_c'].append(np.mean(dict_sample['c_map']))
dict_user['L_m_mass'].append(np.mean(dict_sample['eta_1_map'])+np.mean(dict_sample['eta_2_map'])+np.mean(dict_sample['c_map']))

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
fig.savefig('plot/IC_map_etas_solute.png')
plt.close(fig)

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# PFDEM iteration

dict_sample['Map_known'] = False
dict_sample['i_DEMPF_ite'] = 0
while dict_sample['i_DEMPF_ite'] < dict_user['n_DEMPF_ite']:
    dict_sample['i_DEMPF_ite'] = dict_sample['i_DEMPF_ite'] + 1
    print('\nStep',dict_sample['i_DEMPF_ite'],'/',dict_user['n_DEMPF_ite'],'\n')

    # DEM
    run_yade(dict_user, dict_sample)

    # DEM->PF, PF, PF->DEM
    run_moose(dict_user, dict_sample)

    # compute tracker
    dict_user['L_sum_eta_1'].append(np.sum(dict_sample['eta_1_map']))
    dict_user['L_sum_eta_2'].append(np.sum(dict_sample['eta_2_map']))
    dict_user['L_sum_c'].append(np.sum(dict_sample['c_map']))
    dict_user['L_sum_mass'].append(np.sum(dict_sample['eta_1_map'])+np.sum(dict_sample['eta_2_map'])+np.sum(dict_sample['c_map']))
    dict_user['L_m_eta_1'].append(np.mean(dict_sample['eta_1_map']))
    dict_user['L_m_eta_2'].append(np.mean(dict_sample['eta_2_map']))
    dict_user['L_m_c'].append(np.mean(dict_sample['c_map']))
    dict_user['L_m_mass'].append(np.mean(dict_sample['eta_1_map'])+np.mean(dict_sample['eta_2_map'])+np.mean(dict_sample['c_map']))

    # plot sum eta_i, c
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))
    ax1.plot(dict_user['L_sum_eta_1'])
    ax1.set_title(r'$\Sigma\eta_1$')
    ax2.plot(dict_user['L_sum_eta_2'])
    ax2.set_title(r'$\Sigma\eta_2$')
    ax3.plot(dict_user['L_sum_c'])
    ax3.set_title(r'$\Sigma C$')
    ax4.plot(dict_user['L_sum_mass'])
    ax4.set_title(r'$\Sigma\eta_1 + \Sigma\eta_2 + \Sigma c$')
    fig.savefig('plot/sum_etai_c.png')
    plt.close(fig)

    # plot mean eta_i, c
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))
    ax1.plot(dict_user['L_m_eta_1'])
    ax1.set_title(r'Mean $\eta_1$')
    ax2.plot(dict_user['L_m_eta_2'])
    ax2.set_title(r'Mean $\eta_2$')
    ax3.plot(dict_user['L_m_c'])
    ax3.set_title(r'Mean $c$')
    ax4.plot(dict_user['L_m_mass'])
    ax4.set_title(r'Mean $\eta_1$ + Mean $\eta_2$ + Mean $c$')
    fig.savefig('plot/mean_etai_c.png')
    plt.close(fig)

    # plot performances
    fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
    ax1.plot(dict_user['L_t_dem'], label='DEM')
    ax1.plot(dict_user['L_t_pf'], label='PF')
    ax1.plot(dict_user['L_t_dem_to_pf'], label='DEM to PF')
    ax1.plot(dict_user['L_t_pf_to_dem_1'], label='PF to DEM 1')
    ax1.plot(dict_user['L_t_pf_to_dem_2'], label='PF to DEM 2')
    ax1.legend()
    ax1.set_title('Performances (s)')
    ax1.set_xlabel('Iterations (-)')
    fig.savefig('plot/performances.png')
    plt.close(fig)

    # pp displacement
    L_disp_init = [0]
    L_disp = [0]
    L_strain = [0]
    for i_disp in range(len(dict_user['L_displacement'])):
        L_disp_init.append(L_disp_init[-1]+dict_user['L_displacement'][i_disp])
        if i_disp >= 1:
            L_disp.append(L_disp[-1]+dict_user['L_displacement'][i_disp])
            L_strain.append(L_strain[-1]+dict_user['L_displacement'][i_disp]/(4*dict_user['radius']))
    # compute andrade
    L_andrade = []
    L_strain_log = []
    L_t_log = []
    mean_log_k = 0
    if len(L_strain) > 1:
        for i in range(1,len(L_strain)):
            L_strain_log.append(math.log(abs(L_strain[i])))
            L_t_log.append(math.log(i+1))
            mean_log_k = mean_log_k + (L_strain_log[-1] - 1/3*L_t_log[-1])
        mean_log_k = mean_log_k/len(L_strain) # mean k in Andrade creep law
        # compute fitted Andrade creep law
        for i in range(len(L_t_log)):
            L_andrade.append(mean_log_k + 1/3*L_t_log[i])
    # plot
    fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
    # displacement
    ax1.plot(L_disp)
    ax1.set_title('Displacement (m)')
    # strain
    ax2.plot(L_strain)
    ax2.set_title(r'$\epsilon_y$ (-)')
    ax2.set_xlabel('Times (-)')
    # Andrade
    ax3.plot(L_t_log, L_strain_log)
    ax3.plot(L_t_log, L_andrade, color='k', linestyle='dotted')
    ax3.set_title('Andrade creep law')
    ax3.set_ylabel(r'log(|$\epsilon_y$|) (-)')
    ax3.set_xlabel('log(Times) (-)')
    # close
    fig.savefig('plot/disp_strain_andrade.png')
    plt.close(fig)
    # save
    dict_user['L_disp'] = L_disp
    dict_user['L_disp_init'] = L_disp_init
    dict_user['L_strain'] = L_strain
    dict_user['L_andrade'] = L_andrade
    dict_user['mean_log_k'] = mean_log_k

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# close simulation

# save
with open('dict/dict_user', 'wb') as handle:
    pickle.dump(dict_user, handle, protocol=pickle.HIGHEST_PROTOCOL)
with open('dict/dict_sample', 'wb') as handle:
    pickle.dump(dict_sample, handle, protocol=pickle.HIGHEST_PROTOCOL)

# compute performances
tac = time.perf_counter()
hours = (tac-tic)//(60*60)
minutes = (tac-tic - hours*60*60)//(60)
seconds = int(tac-tic - hours*60*60 - minutes*60)
print("\nSimulation time : "+str(hours)+" hours "+str(minutes)+" minutes "+str(seconds)+" seconds")

print('Simulation ends')

# -*- encoding=utf-8 -*-

from pathlib import Path
import shutil, os, pickle
import numpy as np
import matplotlib.pyplot as plt

# own
from pf_to_dem import *

#------------------------------------------------------------------------------------------------------------------------------------------ #

def create_folder(name):
    '''
    Create a new folder. If it already exists, it is erased.
    '''
    if Path(name).exists():
        shutil.rmtree(name)
    os.mkdir(name)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def tuplet_to_list(tuplet):
    '''
    Convert a tuplet into lists.
    '''
    L_x = []
    L_y = []
    p_center = np.array([0,0])
    n_mean = 0
    for v in tuplet:
        L_x.append(v[0])
        L_y.append(v[1])
        p_center = p_center + np.array([v[0], v[1]])
        n_mean = n_mean + 1
    L_x.append(L_x[0])
    L_y.append(L_y[0])
    p_center = p_center/n_mean
    # translate center to the point (0,0)
    for i in range(len(L_x)):
        L_x[i] = L_x[i] - p_center[0]
        L_y[i] = L_y[i] - p_center[1]
    return L_x, L_y

#------------------------------------------------------------------------------------------------------------------------------------------ #

def tuplet_to_list_no_centerized(tuplet):
    '''
    Convert a tuplet into lists.
    '''
    L_x = []
    L_y = []
    p_center = np.array([0,0])
    n_mean = 0
    for v in tuplet:
        L_x.append(v[0])
        L_y.append(v[1])
        n_mean = n_mean + 1
    L_x.append(L_x[0])
    L_y.append(L_y[0])
    return L_x, L_y

#------------------------------------------------------------------------------------------------------------------------------------------ #

def reduce_n_vtk_files(dict_user, dict_sample):
    '''
    Reduce the number of vtk files for phase-field and dem.

    Warning ! The pf and dem files are not synchronized...
    '''
    if dict_user['n_max_vtk_files'] != None:
        # Phase Field files

        # compute the frequency
        if dict_user['j_total']-1 > dict_user['n_max_vtk_files']:
            f_save = (dict_user['j_total']-1)/(dict_user['n_max_vtk_files']-1)
        else :
            f_save = 1
        # post proccess index
        i_save = 0

        # iterate on time 
        for iteration in range(dict_user['j_total']):
            iteration_str = index_to_str(iteration) # from pf_to_dem.py 
            if iteration >= f_save*i_save:
                i_save_str = index_to_str(i_save) # from pf_to_dem.py
                # rename .pvtu
                os.rename('vtk/pf_'+iteration_str+'.pvtu','vtk/pf_'+i_save_str+'.pvtu')
                # write .pvtu to save all vtk
                file = open('vtk/pf_'+i_save_str+'.pvtu','w')
                file.write('''<?xml version="1.0"?>
                <VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">
                \t<PUnstructuredGrid GhostLevel="1">
                \t\t<PPointData>
                \t\t\t<PDataArray type="Float64" Name="as"/>
                \t\t\t<PDataArray type="Float64" Name="kc"/>
                \t\t\t<PDataArray type="Float64" Name="eta1"/>
                \t\t\t<PDataArray type="Float64" Name="eta2"/>
                \t\t\t<PDataArray type="Float64" Name="c"/>
                \t\t</PPointData>
                \t\t<PCellData>
                \t\t\t<PDataArray type="Int32" Name="libmesh_elem_id"/>
                \t\t\t<PDataArray type="Int32" Name="subdomain_id"/>
                \t\t\t<PDataArray type="Int32" Name="processor_id"/>
                \t\t</PCellData>
                \t\t<PPoints>
                \t\t\t<PDataArray type="Float64" Name="Points" NumberOfComponents="3"/>
                \t\t</PPoints>''')
                line = ''
                for i_proc in range(dict_user['n_proc']):
                    line = line + '''\t\t<Piece Source="pf_'''+i_save_str+'''_'''+str(i_proc)+'''.vtu"/>\n'''
                file.write(line)
                file.write('''\t</PUnstructuredGrid>
                </VTKFile>''')
                file.close()
                # rename .vtk
                for i_proc in range(dict_user['n_proc']):
                    os.rename('vtk/pf_'+iteration_str+'_'+str(i_proc)+'.vtu','vtk/pf_'+i_save_str+'_'+str(i_proc)+'.vtu')
                i_save = i_save + 1 
            else:
                # delete files
                os.remove('vtk/pf_'+iteration_str+'.pvtu')
                for i_proc in range(dict_user['n_proc']):
                    os.remove('vtk/pf_'+iteration_str+'_'+str(i_proc)+'.vtu')
        # .e file
        os.remove('vtk/pf_out.e')
        # other files
        j = 0
        j_str = index_to_str(j)
        filepath = Path('vtk/pf_other_'+j_str+'.pvtu')
        while filepath.exists():
            for i_proc in range(dict_user['n_proc']):
                os.remove('vtk/pf_other_'+j_str+'_'+str(i_proc)+'.vtu')
            os.remove('vtk/pf_other_'+j_str+'.pvtu')
            j = j + 1
            j_str = index_to_str(j)
            filepath = Path('vtk/pf_other_'+j_str+'.pvtu')

        # DEM files

        # compute the frequency
        if 2*dict_user['n_DEMPF_ite']-1 > dict_user['n_max_vtk_files']:
            f_save = (2*dict_user['n_DEMPF_ite']-1)/(dict_user['n_max_vtk_files']-1)
        else :
            f_save = 1
        # post proccess index
        i_save = 0

        # iterate on time 
        for iteration in range(2*dict_user['n_DEMPF_ite']):
            iteration_str = str(iteration) # from pf_to_dem.py 
            if iteration >= f_save*i_save:
                i_save_str = str(i_save) # from pf_to_dem.py
                os.rename('vtk/2grains_'+iteration_str+'.vtk', 'vtk/2grains_'+i_save_str+'.vtk')
                os.rename('vtk/grain1_'+iteration_str+'.vtk', 'vtk/grain1_'+i_save_str+'.vtk')
                os.rename('vtk/grain2_'+iteration_str+'.vtk', 'vtk/grain2_'+i_save_str+'.vtk')
                i_save = i_save + 1
            else :
                os.remove('vtk/2grains_'+iteration_str+'.vtk')
                os.remove('vtk/grain1_'+iteration_str+'.vtk')
                os.remove('vtk/grain2_'+iteration_str+'.vtk')

#------------------------------------------------------------------------------------------------------------------------------------------ #

def save_mesh_database(dict_user, dict_sample):
    '''
    Save mesh database.
    '''
    # creating a database
    if not Path('mesh_map.database').exists():
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'x_min': dict_user['x_min'],
        'x_max': dict_user['x_max'],
        'y_min': dict_user['y_min'],
        'y_max': dict_user['y_max'],
        'n_mesh_x': dict_user['n_mesh_x'],
        'n_mesh_y': dict_user['n_mesh_y'],
        'L_L_i_XYZ_used': dict_sample['L_L_i_XYZ_used'],
        'L_XYZ': dict_sample['L_XYZ']
        }
        dict_database = {'Run_1': dict_data}
        with open('mesh_map.database', 'wb') as handle:
                pickle.dump(dict_database, handle, protocol=pickle.HIGHEST_PROTOCOL)
    # updating a database
    else :
        with open('mesh_map.database', 'rb') as handle:
            dict_database = pickle.load(handle)
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'x_min': dict_user['x_min'],
        'x_max': dict_user['x_max'],
        'y_min': dict_user['y_min'],
        'y_max': dict_user['y_max'],
        'n_mesh_x': dict_user['n_mesh_x'],
        'n_mesh_y': dict_user['n_mesh_y'],
        'L_L_i_XYZ_used': dict_sample['L_L_i_XYZ_used'],
        'L_XYZ': dict_sample['L_XYZ']
        }   
        mesh_map_known = False
        for i_run in range(1,len(dict_database.keys())+1):
            if dict_database['Run_'+str(int(i_run))] == dict_data:
                mesh_map_known = True
        # new entry
        if not mesh_map_known: 
            key_entry = 'Run_'+str(int(len(dict_database.keys())+1))
            dict_database[key_entry] = dict_data
            with open('mesh_map.database', 'wb') as handle:
                pickle.dump(dict_database, handle, protocol=pickle.HIGHEST_PROTOCOL)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def check_mesh_database(dict_user, dict_sample):
    '''
    Check mesh database.
    '''
    if Path('mesh_map.database').exists():
        with open('mesh_map.database', 'rb') as handle:
            dict_database = pickle.load(handle)
        dict_data = {
        'n_proc': dict_user['n_proc'],
        'x_min': dict_user['x_min'],
        'x_max': dict_user['x_max'],
        'y_min': dict_user['y_min'],
        'y_max': dict_user['y_max'],
        'n_mesh_x': dict_user['n_mesh_x'],
        'n_mesh_y': dict_user['n_mesh_y']
        }   
        mesh_map_known = False
        for i_run in range(1,len(dict_database.keys())+1):
            if dict_database['Run_'+str(int(i_run))]['n_proc'] == dict_user['n_proc'] and\
            dict_database['Run_'+str(int(i_run))]['x_min'] == dict_user['x_min'] and\
            dict_database['Run_'+str(int(i_run))]['x_max'] == dict_user['x_max'] and\
            dict_database['Run_'+str(int(i_run))]['y_min'] == dict_user['y_min'] and\
            dict_database['Run_'+str(int(i_run))]['y_max'] == dict_user['y_max'] and\
            dict_database['Run_'+str(int(i_run))]['n_mesh_x'] == dict_user['n_mesh_x'] and\
            dict_database['Run_'+str(int(i_run))]['n_mesh_y'] == dict_user['n_mesh_y'] :
                mesh_map_known = True
                i_known = i_run
        if mesh_map_known :
            dict_sample['Map_known'] = True
            dict_sample['L_L_i_XYZ_used'] = dict_database['Run_'+str(int(i_known))]['L_L_i_XYZ_used']
            dict_sample['L_XYZ'] = dict_database['Run_'+str(int(i_known))]['L_XYZ']
        else :
            dict_sample['Map_known'] = False
    else :
        dict_sample['Map_known'] = False

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_contact_v_s_d(dict_user, dict_sample):
    '''
    Plot figure illustrating the evolution of the volume, surface and height of the contact.
    '''
    # load data
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    # tracker
    dict_user['L_distance_extrema'].append(dict_sample['box_contact_y_max'] - dict_sample['box_contact_y_min'])
    dict_user['L_equivalent_area'].append(1*(dict_sample['box_contact_x_max'] - dict_sample['box_contact_x_min']))
    dict_user['L_contact_volume_box'].append(1*(dict_sample['box_contact_x_max'] - dict_sample['box_contact_x_min'])*(dict_sample['box_contact_y_max'] - dict_sample['box_contact_y_min']))
    dict_user['L_contact_overlap'].append(dict_save['contact_overlap']) # computed by Yade
    dict_user['L_contact_area'].append(dict_save['contact_area']) # computed by Yade
    dict_user['L_contact_volume_yade'].append(dict_save['contact_volume'])
    # plot
    fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(16,9))
    # overlap
    ax1.plot(dict_user['L_distance_extrema'], label='Distance vertices')
    ax1.set_title(r'Contact Box height ($m$)',fontsize=20)
    # surface
    ax2.plot(dict_user['L_equivalent_area'])
    ax2.set_title(r'Contact Box width/surface ($m^2$)',fontsize=20)
    # volume
    ax3.plot(dict_user['L_contact_volume_yade'], label='Yade')
    ax3.plot(dict_user['L_contact_volume_box'], label='Box')
    ax3.legend()
    ax3.set_title(r'Volume ($m^3$)',fontsize=20)
    # close
    plt.suptitle('Contact', fontsize=20)
    fig.savefig('plot/contact_h_s_v.png')
    plt.close(fig)            

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_shape_evolution(dict_user, dict_sample):
    '''
    Plot figure illustrating the evolution of grain shapes.
    '''
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
        fig.tight_layout()
        if dict_user['print_all_shape_evolution']:
            fig.savefig('plot/shape_evolution/'+str(dict_sample['i_DEMPF_ite'])+'.png')
        else:
            fig.savefig('plot/shape_evolution.png')
        plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_n_vertices(dict_user, dict_sample):
    '''
    Plot figure illustrating the number of vertices used in Yade.
    '''
    # load data
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    # tracker
    dict_user['L_n_v_1'].append(dict_save['n_v_1']/2)
    dict_user['L_n_v_2'].append(dict_save['n_v_2']/2)
    dict_user['L_n_v_1_target'].append(dict_save['n_v_1_target']/2)
    dict_user['L_n_v_2_target'].append(dict_save['n_v_2_target']/2)

    # plot
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    ax1.plot(dict_user['L_n_v_1'], label='N vertices g1', color='r')
    ax1.plot(dict_user['L_n_v_1_target'], label='N vertices g1 targetted', color='r', linestyle='dotted')
    ax1.plot(dict_user['L_n_v_2'], label='N vertices g2', color='b')
    ax1.plot(dict_user['L_n_v_2_target'], label='N vertices g2 targetted', color='b', linestyle='dotted')
    ax1.legend()
    fig.tight_layout()
    fig.savefig('plot/n_vertices.png')
    plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_sum_mean_etai_c(dict_user, dict_sample):
    '''
    Plot figure illustrating the sum and the mean of etai and c.
    '''
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
    fig.tight_layout()
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
    fig.tight_layout()
    fig.savefig('plot/mean_etai_c.png')
    plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_performances(dict_user, dict_sample):
    '''
    Plot figure illustrating the time performances of the algorithm.
    '''
    fig, (ax1) = plt.subplots(nrows=1,ncols=1,figsize=(16,9))
    ax1.plot(dict_user['L_t_dem'], label='DEM')
    ax1.plot(dict_user['L_t_pf'], label='PF')
    ax1.plot(dict_user['L_t_dem_to_pf'], label='DEM to PF')
    ax1.plot(dict_user['L_t_pf_to_dem_1'], label='PF to DEM 1')
    ax1.plot(dict_user['L_t_pf_to_dem_2'], label='PF to DEM 2')
    ax1.legend()
    ax1.set_title('Performances (s)')
    ax1.set_xlabel('Iterations (-)')
    fig.tight_layout()
    fig.savefig('plot/performances.png')
    plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_disp_strain_andrade(dict_user, dict_sample):
    '''
    Plot figure illustrating the displacement, the strain and the fit with the Andrade law.
    '''
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
    fig.tight_layout()
    fig.savefig('plot/disp_strain_andrade.png')
    plt.close(fig)
    # save
    dict_user['L_disp'] = L_disp
    dict_user['L_disp_init'] = L_disp_init
    dict_user['L_strain'] = L_strain
    dict_user['L_andrade'] = L_andrade
    dict_user['mean_log_k'] = mean_log_k

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_maps_configuration(dict_user, dict_sample):
    '''
    Plot figure illustrating the current maps of etai and c.
    '''
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
    fig.tight_layout()
    if dict_user['print_all_map_config']:
        fig.savefig('plot/map_etas_solute/'+str(dict_sample['i_DEMPF_ite'])+'.png')
    else:
        fig.savefig('plot/map_etas_solute.png')
    plt.close(fig)


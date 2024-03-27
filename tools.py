# -*- encoding=utf-8 -*-

from pathlib import Path
import shutil, os, pickle, math
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
        'x_min': min(dict_sample['x_L']),
        'x_max': max(dict_sample['x_L']),
        'y_min': min(dict_sample['y_L']),
        'y_max': max(dict_sample['y_L']),
        'n_mesh_x': len(dict_sample['x_L']),
        'n_mesh_y': len(dict_sample['y_L']),
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
        'x_min': min(dict_sample['x_L']),
        'x_max': max(dict_sample['x_L']),
        'y_min': min(dict_sample['y_L']),
        'y_max': max(dict_sample['y_L']),
        'n_mesh_x': len(dict_sample['x_L']),
        'n_mesh_y': len(dict_sample['y_L']),
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
        'x_min': min(dict_sample['x_L']),
        'x_max': max(dict_sample['x_L']),
        'y_min': min(dict_sample['y_L']),
        'y_max': max(dict_sample['y_L']),
        'n_mesh_x': len(dict_sample['x_L']),
        'n_mesh_y': len(dict_sample['y_L'])
        }   
        mesh_map_known = False
        for i_run in range(1,len(dict_database.keys())+1):
            if dict_database['Run_'+str(int(i_run))]['n_proc'] == dict_user['n_proc'] and\
            dict_database['Run_'+str(int(i_run))]['x_min'] == min(dict_sample['x_L']) and\
            dict_database['Run_'+str(int(i_run))]['x_max'] == max(dict_sample['x_L']) and\
            dict_database['Run_'+str(int(i_run))]['y_min'] == min(dict_sample['y_L']) and\
            dict_database['Run_'+str(int(i_run))]['y_max'] == max(dict_sample['y_L']) and\
            dict_database['Run_'+str(int(i_run))]['n_mesh_x'] == len(dict_sample['x_L']) and\
            dict_database['Run_'+str(int(i_run))]['n_mesh_y'] == len(dict_sample['y_L']) :
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
    if 'contact_h_s_v' in dict_user['L_figures']:
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
        if 'shape_evolution' in dict_user['L_figures']:
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
    if 'n_vertices' in dict_user['L_figures']:
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
    if 'sum_etai_c' in dict_user['L_figures']:
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
    if 'mean_etai_c' in dict_user['L_figures']:
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

def compute_mass(dict_user, dict_sample):
    '''
    Compute the mass at a certain time.
     
    Mass is sum of etai and c.
    '''
    # sum of masses
    dict_user['sum_eta_1_tempo'] = np.sum(dict_sample['eta_1_map'])
    dict_user['sum_eta_2_tempo'] = np.sum(dict_sample['eta_2_map'])
    dict_user['sum_c_tempo'] = np.sum(dict_sample['c_map'])
    dict_user['sum_mass_tempo'] = np.sum(dict_sample['eta_1_map'])+np.sum(dict_sample['eta_2_map'])+np.sum(dict_sample['c_map'])
    
#------------------------------------------------------------------------------------------------------------------------------------------ #

def compute_mass_loss(dict_user, dict_sample, tracker_key):
    '''
    Compute the mass loss from the previous compute_mass() call.
     
    Plot in the given tracker.
    Mass is sum of etai and c.
    '''
    # delta masses
    deta1 = np.sum(dict_sample['eta_1_map']) - dict_user['sum_eta_1_tempo']
    deta2 = np.sum(dict_sample['eta_2_map']) - dict_user['sum_eta_2_tempo']
    dc = np.sum(dict_sample['c_map']) - dict_user['sum_c_tempo']
    dm = np.sum(dict_sample['eta_1_map'])+np.sum(dict_sample['eta_2_map'])+np.sum(dict_sample['c_map']) - dict_user['sum_mass_tempo']
    
    # save
    dict_user[tracker_key+'_eta1'].append(deta1)
    dict_user[tracker_key+'_eta2'].append(deta2)
    dict_user[tracker_key+'_c'].append(dc)
    dict_user[tracker_key+'_m'].append(dm)

    # plot
    if 'mass_loss' in dict_user['L_figures']:
        fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2,ncols=2,figsize=(16,9))
        ax1.plot(dict_user[tracker_key+'_eta1'])
        ax1.set_title(r'$\eta_1$ loss')
        ax2.plot(dict_user[tracker_key+'_eta2'])
        ax2.set_title(r'$\eta_2$ loss')
        ax3.plot(dict_user[tracker_key+'_c'])
        ax3.set_title(r'$c$ loss')
        ax4.plot(dict_user[tracker_key+'_m'])
        ax4.set_title(r'$\eta_1$ + $\eta_2$ + $c$ loss')
        fig.tight_layout()
        fig.savefig('plot/'+tracker_key+'.png')
        plt.close(fig)

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_performances(dict_user, dict_sample):
    '''
    Plot figure illustrating the time performances of the algorithm.
    '''
    if 'performances' in dict_user['L_figures']:
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
    # pp time PF
    for i_dt in range(len(dict_user['L_dt_PF'])):
        if i_dt == 0:
            L_t_PF = [dict_user['L_dt_PF'][0]]
        else:
            L_t_PF.append(L_t_PF[-1] + dict_user['L_dt_PF'][i_dt])

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
            L_t_log.append(math.log(L_t_PF[i-1]))
            mean_log_k = mean_log_k + (L_strain_log[-1] - 1/3*L_t_log[-1])
        mean_log_k = mean_log_k/len(L_strain) # mean k in Andrade creep law
        # compute fitted Andrade creep law
        for i in range(len(L_t_log)):
            L_andrade.append(mean_log_k + 1/3*L_t_log[i])
    # plot
    if 'disp_strain_andrade' in dict_user['L_figures'] and dict_sample['i_DEMPF_ite'] > 10:
        fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
        # displacement
        ax1.plot(L_t_PF, L_disp)
        ax1.set_title('Displacement (m)')
        ax1.set_xlabel('PF Times (-)')
        # strain
        ax2.plot(L_t_PF, L_strain)
        ax2.set_title(r'$\epsilon_y$ (-)')
        ax2.set_xlabel('PF Times (-)')
        # Andrade
        ax3.plot(L_t_log, L_strain_log)
        ax3.plot(L_t_log, L_andrade, color='k', linestyle='dotted')
        ax3.set_title('Andrade creep law')
        ax3.set_ylabel(r'log(|$\epsilon_y$|) (-)')
        ax3.set_xlabel('log(PF Times) (-)')
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

def plot_sample_height(dict_user, dict_sample):
    '''
    Plot figure illustrating the displacement, the strain and the fit with the Andrade law, assuming the sample height.
    '''
    # pp time PF
    for i_dt in range(len(dict_user['L_dt_PF'])):
        if i_dt == 0:
            L_t_PF = [dict_user['L_dt_PF'][0]]
        else:
            L_t_PF.append(L_t_PF[-1] + dict_user['L_dt_PF'][i_dt])

    # pp sample height
    L_strain = []
    for i_height in range(len(dict_user['L_sample_height'])):
        L_strain.append((dict_user['L_sample_height'][i_height]-4*dict_user['radius'])/(4*dict_user['radius']))

    # compute andrade
    L_andrade = []
    L_strain_log = []
    L_t_log = []
    mean_log_k = 0
    if len(L_strain) > 1:
        for i in range(1,len(L_strain)):
            L_strain_log.append(math.log(abs(L_strain[i])))
            L_t_log.append(math.log(L_t_PF[i-1]))
            mean_log_k = mean_log_k + (L_strain_log[-1] - 1/3*L_t_log[-1])
        mean_log_k = mean_log_k/len(L_strain) # mean k in Andrade creep law
        # compute fitted Andrade creep law
        for i in range(len(L_t_log)):
            L_andrade.append(mean_log_k + 1/3*L_t_log[i])
    # plot
    if 'sample_height' in dict_user['L_figures'] and dict_sample['i_DEMPF_ite'] > 10:
        fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
        # displacement
        ax1.plot(L_t_PF[1:], dict_user['L_sample_height'][1:])
        ax1.set_title('Sample height (m)')
        ax1.set_xlabel('PF Times (-)')
        # strain
        ax2.plot(L_t_PF, L_strain)
        ax2.set_title(r'$\epsilon_y$ (-)')
        ax2.set_xlabel('PF Times (-)')
        # Andrade
        ax3.plot(L_t_log, L_strain_log)
        ax3.plot(L_t_log, L_andrade, color='k', linestyle='dotted')
        ax3.set_title('Andrade creep law')
        ax3.set_ylabel(r'log(|$\epsilon_y$|) (-)')
        ax3.set_xlabel('log(Times) (-)')
        # close
        fig.tight_layout()
        fig.savefig('plot/sample_height.png')
        plt.close(fig)
    # save
    dict_user['L_strain'] = L_strain
    dict_user['L_andrade'] = L_andrade
    dict_user['mean_log_k'] = mean_log_k

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_y_contactPoint(dict_user, dict_sample):
    '''
    Plot figure illustrating the displacement, the strain and the fit with the Andrade law, assuming the y coordinate of the contact point.
    '''
    # pp time PF
    for i_dt in range(len(dict_user['L_dt_PF'])):
        if i_dt == 0:
            L_t_PF = [dict_user['L_dt_PF'][0]]
        else:
            L_t_PF.append(L_t_PF[-1] + dict_user['L_dt_PF'][i_dt])

    # pp sample height
    L_strain = []
    for i_height in range(0,len(dict_user['L_y_contactPoint'])):
        L_strain.append((dict_user['L_y_contactPoint'][i_height]-dict_user['L_y_contactPoint'][0])/(4*dict_user['radius']))

    # compute andrade
    L_andrade = []
    L_strain_log = []
    L_t_log = []
    mean_log_k = 0
    if len(L_strain) > 1:
        for i in range(1,len(L_strain)):
            L_strain_log.append(math.log(abs(L_strain[i])))
            L_t_log.append(math.log(L_t_PF[i-1]))
            mean_log_k = mean_log_k + (L_strain_log[-1] - 1/3*L_t_log[-1])
        mean_log_k = mean_log_k/len(L_strain) # mean k in Andrade creep law
        # compute fitted Andrade creep law
        for i in range(len(L_t_log)):
            L_andrade.append(mean_log_k + 1/3*L_t_log[i])
    # plot
    if 'y_contactPoint' in dict_user['L_figures'] and dict_sample['i_DEMPF_ite'] > 10:
        fig, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3,figsize=(16,9))
        # displacement
        ax1.plot(L_t_PF, dict_user['L_y_contactPoint'])
        ax1.set_title('Y coordinate of the contact point (m)')
        ax1.set_xlabel('PF Times (-)')
        # strain
        ax2.plot(L_t_PF, L_strain)
        ax2.set_title(r'$\epsilon_y$ (-)')
        ax2.set_xlabel('PF Times (-)')
        # Andrade
        ax3.plot(L_t_log, L_strain_log)
        ax3.plot(L_t_log, L_andrade, color='k', linestyle='dotted')
        ax3.set_title('Andrade creep law')
        ax3.set_ylabel(r'log(|$\epsilon_y$|) (-)')
        ax3.set_xlabel('log(Times) (-)')
        # close
        fig.tight_layout()
        fig.savefig('plot/y_contact_point.png')
        plt.close(fig)
    # save
    dict_user['L_strain'] = L_strain
    dict_user['L_andrade'] = L_andrade
    dict_user['mean_log_k'] = mean_log_k

#------------------------------------------------------------------------------------------------------------------------------------------ #

def plot_maps_configuration(dict_user, dict_sample):
    '''
    Plot figure illustrating the current maps of etai and c.
    '''
    # Plot
    if 'maps' in dict_user['L_figures']:
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

#------------------------------------------------------------------------------------------------------------------------------------------ #

def compute_sphericities(L_vertices):
    '''
    Compute sphericity of the particle with five parameters.

    The parameters used are the area, the diameter, the circle ratio, the perimeter and the width to length ratio sphericity.
    See Zheng, J., Hryciw, R.D. (2015) Traditional soil particle sphericity, roundness and surface roughness by computational geometry, Geotechnique, Vol 65
    '''
    # adapt list
    L_vertices_x, L_vertices_y = tuplet_to_list_no_centerized(L_vertices)
    L_vertices = []
    for i_v in range(len(L_vertices_x)):
        L_vertices.append(np.array([L_vertices_x[i_v], L_vertices_y[i_v]]))

    #Find the minimum circumscribing circle
    #look for the two farthest and nearest points
    MaxDistance = 0
    for i_p in range(0,len(L_vertices)-2):
        for j_p in range(i_p+1,len(L_vertices)-1):
            Distance = np.linalg.norm(L_vertices[i_p]-L_vertices[j_p])
            if Distance > MaxDistance :
                ij_farthest = (i_p,j_p)
                MaxDistance = Distance

    #Trial circle
    center_circumscribing = (L_vertices[ij_farthest[0]]+L_vertices[ij_farthest[1]])/2
    radius_circumscribing = MaxDistance/2
    Circumscribing_Found = True
    Max_outside_distance = radius_circumscribing
    for i_p in range(len(L_vertices)-1):
        #there is a margin here because of the numerical approximation
        if np.linalg.norm(L_vertices[i_p]-center_circumscribing) > (1+0.05)*radius_circumscribing and i_p not in ij_farthest: #vertex outside the trial circle
            Circumscribing_Found = False
            if np.linalg.norm(L_vertices[i_p]-center_circumscribing) > Max_outside_distance:
                k_outside_farthest = i_p
                Max_outside_distance = np.linalg.norm(L_vertices[i_p]-center_circumscribing)
    #The trial guess does not work
    if not Circumscribing_Found:
        L_ijk_circumscribing = [ij_farthest[0],ij_farthest[1],k_outside_farthest]
        center_circumscribing, radius_circumscribing = FindCircleFromThreePoints(L_vertices[L_ijk_circumscribing[0]],L_vertices[L_ijk_circumscribing[1]],L_vertices[L_ijk_circumscribing[2]])
        Circumscribing_Found = True
        for i_p in range(len(L_vertices)-1):
            #there is 1% margin here because of the numerical approximation
            if np.linalg.norm(L_vertices[i_p]-center_circumscribing) > (1+0.05)*radius_circumscribing and i_p not in L_ijk_circumscribing: #vertex outside the circle computed
                Circumscribing_Found = False

    #look for length and width
    length = MaxDistance
    u_maxDistance = (L_vertices[ij_farthest[0]]-L_vertices[ij_farthest[1]])/np.linalg.norm(L_vertices[ij_farthest[0]]-L_vertices[ij_farthest[1]])
    v_maxDistance = np.array([u_maxDistance[1], -u_maxDistance[0]])
    MaxWidth = 0
    for i_p in range(0,len(L_vertices)-2):
        for j_p in range(i_p+1,len(L_vertices)-1):
            Distance = abs(np.dot(L_vertices[i_p]-L_vertices[j_p],v_maxDistance))
            if Distance > MaxWidth :
                ij_width = (i_p,j_p)
                MaxWidth = Distance
    width = MaxWidth

    #look for maximum inscribed circle
    #discretisation of the grain
    l_x_inscribing = np.linspace(min(L_vertices_x),max(L_vertices_x), 100)
    l_y_inscribing = np.linspace(min(L_vertices_y),max(L_vertices_y), 100)
    #creation of an Euclidean distance map to the nearest boundary vertex
    map_inscribing = np.zeros((100, 100))
    #compute the map
    for i_x in range(100):
        for i_y in range(100):
            p = np.array([l_x_inscribing[i_x], l_y_inscribing[-1-i_y]])
            #work only if the point is inside the grain
            if P_is_inside(L_vertices, p):
                #look for the nearest vertex
                MinDistance = None
                for q in L_vertices[:-1]:
                    Distance = np.linalg.norm(p-q)
                    if MinDistance == None or Distance < MinDistance:
                        MinDistance = Distance
                map_inscribing[-1-i_y, i_x] = MinDistance
            else :
                map_inscribing[-1-i_y, i_x] = 0
    #look for the peak of the map
    index_max = np.argmax(map_inscribing)
    l = index_max//100
    c = index_max%100
    radius_inscribing = map_inscribing[l, c]

    #Compute surface of the grain 
    #Sinus law
    meanPoint = np.mean(L_vertices[:-1], axis=0)
    SurfaceParticle = 0
    for i_triangle in range(len(L_vertices)-1):
        AB = np.array(L_vertices[i_triangle]-meanPoint)
        AC = np.array(L_vertices[i_triangle+1]-meanPoint)
        SurfaceParticle = SurfaceParticle + 0.5*np.linalg.norm(np.cross(AB, AC))

    #Area Sphericity
    if Circumscribing_Found :
        SurfaceCircumscribing = math.pi*radius_circumscribing**2
        AreaSphericity = SurfaceParticle / SurfaceCircumscribing
    else :
        AreaSphericity = 1

    #Diameter Sphericity
    if Circumscribing_Found :
        DiameterSameAreaParticle = 2*math.sqrt(SurfaceParticle/math.pi)
        DiameterCircumscribing = radius_circumscribing*2
        DiameterSphericity = DiameterSameAreaParticle / DiameterCircumscribing
    else :
        DiameterSphericity = 1

    #Circle Ratio Sphericity
    if Circumscribing_Found :
        DiameterInscribing = radius_inscribing*2
        CircleRatioSphericity = DiameterInscribing / DiameterCircumscribing
    else : 
        CircleRatioSphericity = 1
    
    #Perimeter Sphericity
    PerimeterSameAreaParticle = 2*math.sqrt(SurfaceParticle*math.pi)
    PerimeterParticle = 0
    for i in range(len(L_vertices)-1):
        PerimeterParticle = PerimeterParticle + np.linalg.norm(L_vertices[i+1]-L_vertices[i])
    PerimeterSphericity = PerimeterSameAreaParticle / PerimeterParticle

    #Width to length ratio Spericity
    WidthToLengthRatioSpericity = width / length

    return AreaSphericity, DiameterSphericity, CircleRatioSphericity, PerimeterSphericity, WidthToLengthRatioSpericity

#------------------------------------------------------------------------------------------------------------------------------------------ #

def P_is_inside(L_vertices, P):
    '''
    Determine if a point P is inside of a grain

    Make a slide on constant y. Every time a border is crossed, the point switches between in and out.
    see Franklin 1994, see Alonso-Marroquin 2009
    '''
    counter = 0
    for i_p_border in range(len(L_vertices)-1):
        #consider only points if the coordinates frame the y-coordinate of the point
        if (L_vertices[i_p_border][1]-P[1])*(L_vertices[i_p_border+1][1]-P[1]) < 0 :
            x_border = L_vertices[i_p_border][0] + (L_vertices[i_p_border+1][0]-L_vertices[i_p_border][0])*(P[1]-L_vertices[i_p_border][1])/(L_vertices[i_p_border+1][1]-L_vertices[i_p_border][1])
            if x_border > P[0] :
                counter = counter + 1
    if counter % 2 == 0:
        return False
    else :
        return True
    
#------------------------------------------------------------------------------------------------------------------------------------------ #

def FindCircleFromThreePoints(P1, P2, P3):
    '''
    Compute the circumscribing circle of a triangle defined by three points.

    https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/
    '''
    # Line P1P2 is represented as ax + by = c and line P2P3 is represented as ex + fy = g
    a, b, c = lineFromPoints(P1, P2)
    e, f, g = lineFromPoints(P2, P3)

    # Converting lines P1P2 and P2P3 to perpendicular bisectors.
    #After this, L : ax + by = c and M : ex + fy = g
    a, b, c = perpendicularBisectorFromLine(P1, P2, a, b, c)
    e, f, g = perpendicularBisectorFromLine(P2, P3, e, f, g)

    # The point of intersection of L and M gives the circumcenter
    circumcenter = lineLineIntersection(a, b, c, e, f, g)

    if np.linalg.norm(circumcenter - np.array([10**9,10**9])) == 0:
        raise ValueError('The given points do not form a triangle and are collinear...')
    else :
        #compute the radius
        radius = max([np.linalg.norm(P1-circumcenter), np.linalg.norm(P2-circumcenter), np.linalg.norm(P3-circumcenter)])

    return circumcenter, radius

#------------------------------------------------------------------------------------------------------------------------------------------ #

def lineFromPoints(P, Q):
    '''
    Function to find the line given two points

    Used in FindCircleFromThreePoints().
    The equation is c = ax + by.
    https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/
    '''
    a = Q[1] - P[1]
    b = P[0] - Q[0]
    c = a * (P[0]) + b * (P[1])
    return a, b, c

#------------------------------------------------------------------------------------------------------------------------------------------ #

def lineLineIntersection(a1, b1, c1, a2, b2, c2):
    '''
    Returns the intersection point of two lines.

    Used in FindCircleFromThreePoints().
    https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/
    '''
    determinant = a1 * b2 - a2 * b1
    if (determinant == 0):
        # The lines are parallel.
        return np.array([10**9,10**9])
    else:
        x = (b2 * c1 - b1 * c2)//determinant
        y = (a1 * c2 - a2 * c1)//determinant
        return np.array([x, y])

#------------------------------------------------------------------------------------------------------------------------------------------ #

def perpendicularBisectorFromLine(P, Q, a, b, c):
    '''
    Function which converts the input line to its perpendicular bisector.

    Used in FindCircleFromThreePoints().
    The equation is c = ax + by.
    https://www.geeksforgeeks.org/program-find-circumcenter-triangle-2/
    '''
    mid_point = [(P[0] + Q[0])//2, (P[1] + Q[1])//2]
    # c = -bx + ay
    c = -b * (mid_point[0]) + a * (mid_point[1])
    temp = a
    a = -b
    b = temp
    return a, b, c

#------------------------------------------------------------------------------------------------------------------------------------------ #

def remesh(dict_user, dict_sample):
    '''
    Remesh the problem.
    
    Eta1, Eta2, c maps are updated
    x_L, n_mesh_x, y_L, n_mesh_y are updated.
    '''
    # search the grain boundaries
    y_min = dict_sample['y_L'][-1]
    y_max = dict_sample['y_L'][0]
    x_min = dict_sample['x_L'][-1]
    x_max = dict_sample['x_L'][0]
    # iterate on y
    for i_y in range(len(dict_sample['y_L'])):
        if max(dict_sample['eta_1_map'][-1-i_y, :]) > 0.5 or\
           max(dict_sample['eta_2_map'][-1-i_y, :]) > 0.5:
            if dict_sample['y_L'][i_y] < y_min : 
                y_min = dict_sample['y_L'][i_y]
            if dict_sample['y_L'][i_y] > y_max :
                y_max = dict_sample['y_L'][i_y]
    # iterate on x
    for i_x in range(len(dict_sample['x_L'])):
        if max(dict_sample['eta_1_map'][:, i_x]) > 0.5 or\
           max(dict_sample['eta_2_map'][:, i_x]) > 0.5:
            if dict_sample['x_L'][i_x] < x_min : 
                x_min = dict_sample['x_L'][i_x]
            if dict_sample['x_L'][i_x] > x_max :
                x_max = dict_sample['x_L'][i_x]
    # compute the domain boundaries (grain boundaries + margins)
    x_min_dom = x_min - dict_user['margin_mesh_domain']
    x_max_dom = x_max + dict_user['margin_mesh_domain']
    y_min_dom = y_min - dict_user['margin_mesh_domain']
    y_max_dom = y_max + dict_user['margin_mesh_domain']
    # compute the new x_L and y_L
    x_L = np.arange(x_min_dom, x_max_dom, dict_user['size_x_mesh'])
    n_mesh_x = len(x_L)
    y_L = np.arange(y_min_dom, y_max_dom, dict_user['size_y_mesh'])
    n_mesh_y = len(y_L)
    delta_x_max = 0
    delta_y_max = 0
    # compute the new maps
    eta_1_map = np.zeros((n_mesh_y, n_mesh_x))
    eta_2_map = np.zeros((n_mesh_y, n_mesh_x))
    c_map = np.ones((n_mesh_y, n_mesh_x))
    # iterate on lines
    for i_y in range(len(y_L)):
        # addition
        if y_L[i_y] < dict_sample['y_L'][0] or \
           dict_sample['y_L'][-1] < y_L[i_y]:
            # iterate on columns
            for i_x in range(len(x_L)):
                eta_1_map[-1-i_y, i_x] = 0
                eta_2_map[-1-i_y, i_x] = 0
                c_map[-1-i_y, i_x] = 1
        # extraction
        else :     
            # iterate on columns
            for i_x in range(len(x_L)):
                # addition
                if x_L[i_x] < dict_sample['x_L'][0] or \
                   dict_sample['x_L'][-1] < x_L[i_x]: 
                    eta_1_map[-1-i_y, i_x] = 0
                    eta_2_map[-1-i_y, i_x] = 0
                    c_map[-1-i_y, i_x] = 1
                # extraction
                else :
                    # find nearest node to old node
                    L_search = list(abs(np.array(dict_sample['x_L']-x_L[i_x])))
                    i_x_old = L_search.index(min(L_search))
                    delta_x = min(L_search)
                    L_search = list(abs(np.array(dict_sample['y_L']-y_L[i_y])))
                    i_y_old = L_search.index(min(L_search))
                    delta_y = min(L_search)
                    # track
                    if delta_x > delta_x_max:
                        delta_x_max = delta_x
                    if delta_y > delta_y_max:
                        delta_y_max = delta_y
                    # update
                    eta_1_map[-1-i_y, i_x] = dict_sample['eta_1_map'][-1-i_y_old, i_x_old]  
                    eta_2_map[-1-i_y, i_x] = dict_sample['eta_2_map'][-1-i_y_old, i_x_old]  
                    c_map[-1-i_y, i_x] = dict_sample['c_map'][-1-i_y_old, i_x_old]  
    # tracking
    dict_user['L_x_min_dom'].append(min(x_L))
    dict_user['L_x_max_dom'].append(max(x_L))
    dict_user['L_y_min_dom'].append(min(y_L))
    dict_user['L_y_max_dom'].append(max(y_L))
    dict_user['L_delta_x_max'].append(delta_x_max)
    dict_user['L_delta_y_max'].append(delta_y_max)
    # plot 
    if 'dim_dom' in dict_user['L_figures']:
        fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
        ax1.plot(dict_user['L_x_min_dom'], label='x_min')
        ax1.plot(dict_user['L_x_max_dom'], label='x_max')
        ax1.plot(dict_user['L_y_min_dom'], label='y_min')
        ax1.plot(dict_user['L_y_max_dom'], label='y_max')
        ax1.legend(fontsize=20)
        ax1.set_title(r'Domain Dimensions',fontsize = 30)
        fig.tight_layout()
        fig.savefig('plot/dim_dom.png')
        plt.close(fig)

        fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
        ax1.plot(dict_user['L_delta_x_max'], label=r'max $\Delta$x')
        ax1.plot(dict_user['L_delta_y_max'], label=r'max $\Delta$y')
        ax1.legend(fontsize=20)
        ax1.set_title(r'Interpolation errors',fontsize = 30)
        fig.tight_layout()
        fig.savefig('plot/dim_dom_delta.png')
        plt.close(fig)
    # save 
    dict_sample['x_L'] = x_L
    dict_user['n_mesh_x'] = n_mesh_x
    dict_sample['y_L'] = y_L
    dict_user['n_mesh_y'] = n_mesh_y
    dict_sample['eta_1_map'] = eta_1_map
    dict_sample['eta_2_map'] = eta_2_map
    dict_sample['c_map'] = c_map

                
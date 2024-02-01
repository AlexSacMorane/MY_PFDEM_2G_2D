# -*- encoding=utf-8 -*-

import pickle, math, os, shutil
from pathlib import Path
from scipy.ndimage import binary_dilation
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------#

def move_phasefield(dict_user, dict_sample):
    '''
    Move phase field maps by interpolation.
    '''
    # g2 (as g1 is fixed)

    # load data
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    displacement_g2 = dict_save['displacement_g2']

    # tracker
    dict_user['L_displacement'].append(displacement_g2[1])

    # loading old variables
    eta_1_map = dict_sample['eta_1_map']
    eta_2_map = dict_sample['eta_2_map']
    c_map = dict_sample['c_map']

    # updating phase map
    print('Updating phase field maps')
    eta_1_map_new = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    eta_2_map_new = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    c_map_new = c_map # not updated

    # iteration on y
    i_y_old = 0
    for i_y in range(len(dict_sample['y_L'])):
        y = dict_sample['y_L'][i_y]

        # eta 1, fixed
        eta_1_map_new = eta_1_map

        # eta 2
        if displacement_g2[1] < 0:
            if y-displacement_g2[1] <= dict_sample['y_L'][-1]:
                # look for window
                while not (dict_sample['y_L'][i_y_old] <= y-displacement_g2[1] and y-displacement_g2[1] < dict_sample['y_L'][i_y_old+1]):
                    i_y_old = i_y_old + 1
                # interpolate
                eta_2_map_new[-1-i_y, :] = (eta_2_map[-1-(i_y_old+1), :] - eta_2_map[-1-i_y_old, :])/(dict_sample['y_L'][i_y_old+1] - dict_sample['y_L'][i_y_old])*\
                                           (y-displacement_g2[1] - dict_sample['y_L'][i_y_old]) + eta_2_map[-1-i_y_old, :]
        elif displacement_g2[1] > 0:
            if dict_sample['y_L'][0] <= y-displacement_g2[1]:
                # look for window
                while not (dict_sample['y_L'][i_y_old] <= y-displacement_g2[1] and y-displacement_g2[1] < dict_sample['y_L'][i_y_old+1]):
                    i_y_old = i_y_old + 1
                # interpolate
                eta_2_map_new[-1-i_y, :] = (eta_2_map[-1-(i_y_old+1), :] - eta_2_map[-1-i_y_old, :])/(dict_sample['y_L'][i_y_old+1] - dict_sample['y_L'][i_y_old])*\
                                           (y-displacement_g2[1] - dict_sample['y_L'][i_y_old]) + eta_2_map[-1-i_y_old, :]
        else :
            eta_2_map_new = eta_2_map

        # solute
        # TO DO : push solute with grain displacement

    # update variables
    dict_sample['eta_1_map'] = eta_1_map_new
    dict_sample['eta_2_map'] = eta_2_map_new
    dict_sample['c_map'] = c_map_new

    # write txts for phase field
    write_eta_txt(dict_user, dict_sample) # phase field
    write_c_txt(dict_user, dict_sample) # solute

# -----------------------------------------------------------------------------#

def compute_ed(dict_user, dict_sample):
    '''
    Compute sink/source term.
    '''
    # load data from dem
    with open('data/dem_to_main.data', 'rb') as handle:
        dict_save = pickle.load(handle)
    contact_area = dict_save['equivalent_area']
    normalForce = dict_save['normalForce']
    contactPoint = dict_save['contact_point']

    # tracker
    dict_user['L_P_applied'].append(np.linalg.norm(normalForce)/contact_area)

    # plot
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    # overlap
    ax1.plot(dict_user['L_P_applied'])
    ax1.set_title('Pressure at the contact (Pa)',fontsize=20)
    fig.savefig('plot/contact_pressure.png')
    plt.close(fig)

    # init
    dict_sample['ed_map'] = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    # cst
    R_gas = 82.06e5 # cm3 Pa K-1 mol-1
    Temp = 25+278   # K
    V_m = 27.1      # cm3 mol-1
    conv = 1        # = 1/C_ref, normalize the sink/source
    s_mesh = 1*(dict_sample['x_L'][1]-dict_sample['x_L'][0]) # surface of a mesh, m2
    v_mesh = 1*(dict_sample['x_L'][1]-dict_sample['x_L'][0])*(dict_sample['y_L'][1]-dict_sample['y_L'][0]) # volume of a mesh, m3

    # iterate on mesh
    for i_x in range(len(dict_sample['x_L'])):
        for i_y in range(len(dict_sample['y_L'])):
            # determine pressure
            if dict_sample['eta_1_map'][i_y, i_x] > 0.5 and dict_sample['eta_2_map'][i_y, i_x] > 0.5: # in the contact
                P =  np.linalg.norm(normalForce)/contact_area # Pa
            else : # not in the contact
                P = 0 # Pa
            r = s_mesh*math.exp(P*V_m/(R_gas*Temp))*(1-dict_sample['c_map'][i_y, i_x]/(dict_user['C_eq']*math.exp(P*V_m/(R_gas*Temp))))
            if r >= 0 :
                r = dict_user['k_diss']*r # dissolution rate
            else :
                r = dict_user['k_prec']*r # precipitation rate
            # the rate is in mol s-1
            r = r/v_mesh # convert in mol m-3 s-1
            r = r*conv # convert in C_ref m-3 s-1
            r = r*2/3 # convert in C_ref m-3
            # save in the map
            dict_sample['ed_map'][i_y, i_x] = r

    # tracker
    dict_user['L_m_ed'].append(np.mean(dict_sample['ed_map']))

    # plot
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    # overlap
    ax1.plot(dict_user['L_m_ed'])
    ax1.set_title('Mean tilting factor (-)',fontsize=20)
    fig.savefig('plot/m_ed.png')
    plt.close(fig)

    # find nearest node to contact point
    L_search = list(abs(np.array(dict_sample['x_L']-contactPoint[0])))
    i_x = L_search.index(min(L_search))
    L_search = list(abs(np.array(dict_sample['y_L']-contactPoint[1])))
    i_y = L_search.index(min(L_search))

    # tracker
    dict_user['L_ed_contact_point'].append(dict_sample['ed_map'][-1-i_y, i_x])

    # plot
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    # overlap
    ax1.plot(dict_user['L_ed_contact_point'])
    ax1.set_title('Tilting factor at contact point (-)',fontsize=20)
    fig.savefig('plot/contact_ed.png')
    plt.close(fig)

    # write ed
    write_ed_txt(dict_user, dict_sample)

#-------------------------------------------------------------------------------

def write_eta_txt(dict_user, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    This .txt file represent the phase field maps.
    '''
    file_to_write_1 = open('data/eta_1.txt','w')
    file_to_write_2 = open('data/eta_2.txt','w')
    # x
    file_to_write_1.write('AXIS X\n')
    file_to_write_2.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    file_to_write_2.write(line)
    # y
    file_to_write_1.write('AXIS Y\n')
    file_to_write_2.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    file_to_write_2.write(line)
    # data
    file_to_write_1.write('DATA\n')
    file_to_write_2.write('DATA\n')
    for j in range(len(dict_sample['y_L'])):
        for i in range(len(dict_sample['x_L'])):
            # grain 1
            file_to_write_1.write(str(dict_sample['eta_1_map'][-1-j,i])+'\n')
            # grain 2
            file_to_write_2.write(str(dict_sample['eta_2_map'][-1-j,i])+'\n')
    # close
    file_to_write_1.close()
    file_to_write_2.close()

#-------------------------------------------------------------------------------

def write_c_txt(dict_user, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    This .txt represents the solute map.
    '''
    file_to_write = open('data/c.txt','w')
    # x
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # y
    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # data
    file_to_write.write('DATA\n')
    for j in range(len(dict_sample['y_L'])):
        for i in range(len(dict_sample['x_L'])):
            file_to_write.write(str(dict_sample['c_map'][-1-j,i])+'\n')
    # close
    file_to_write.close()

#-------------------------------------------------------------------------------

def write_ed_txt(dict_user, dict_sample):
    '''
    Write a .txt file needed for MOOSE simulation.

    This .txt represents the map of the external mechanical energy.
    '''
    file_to_write = open('data/ed.txt','w')
    # x
    file_to_write.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # y
    file_to_write.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # data
    file_to_write.write('DATA\n')
    for j in range(len(dict_sample['y_L'])):
        for i in range(len(dict_sample['x_L'])):
            file_to_write.write(str(dict_sample['ed_map'][-1-j,i])+'\n')
    # close
    file_to_write.close()

#-------------------------------------------------------------------------------

def compute_kc(dict_user, dict_sample):
    '''
    Compute the diffusion coefficient of the solute.
    Then Write a .txt file needed for MOOSE simulation.

    This .txt file represent the phase field maps.
    '''
    # compute
    kc_map = np.array(np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x'])), dtype = bool)
    # iterate on x and y
    for i_y in range(len(dict_sample['y_L'])):
        for i_x in range(len(dict_sample['x_L'])):
            if dict_sample['eta_1_map'][i_y, i_x] < 0.5 and dict_sample['eta_2_map'][i_y, i_x] < 0.5: # out of the grain
                kc_map[i_y, i_x] = True
            elif dict_sample['eta_1_map'][i_y, i_x] > 0.5 and dict_sample['eta_2_map'][i_y, i_x] > 0.5: # in the contact
                kc_map[i_y, i_x] = True
            else :
                kc_map[i_y, i_x] = False

    # dilation
    dilated_M = binary_dilation(kc_map, dict_user['struct_element'])

    #compute the map of the solute diffusion coefficient
    kc_map = np.array(np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x'])))
    for i_y in range(len(dict_sample['y_L'])):
        for i_x in range(len(dict_sample['x_L'])):
            if dilated_M[i_y, i_x]:
                kc_map[i_y, i_x] = dict_user['D_solute']

    # write
    file_to_write_1 = open('data/kc.txt','w')
    # x
    file_to_write_1.write('AXIS X\n')
    line = ''
    for x in dict_sample['x_L']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    # y
    file_to_write_1.write('AXIS Y\n')
    line = ''
    for y in dict_sample['y_L']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    # data
    file_to_write_1.write('DATA\n')
    for j in range(len(dict_sample['y_L'])):
        for i in range(len(dict_sample['x_L'])):
            # grain 1
            file_to_write_1.write(str(kc_map[-1-j,i])+'\n')
    # close
    file_to_write_1.close()

#-------------------------------------------------------------------------------

def write_i(dict_user, dict_sample):
  '''
  Create the .i file to run MOOSE simulation.

  The file is generated from a template nammed PF_ACS_base.i
  '''
  file_to_write = open('pf.i','w')
  file_to_read = open('pf_base.i','r')
  lines = file_to_read.readlines()
  file_to_read.close()

  j = 0
  for line in lines :
    j = j + 1
    if j == 4:
      line = line[:-1] + ' ' + str(len(dict_sample['x_L'])-1)+'\n'
    elif j == 5:
      line = line[:-1] + ' ' + str(len(dict_sample['y_L'])-1)+'\n'
    elif j == 6:
      line = line[:-1] + ' ' + str(min(dict_sample['x_L']))+'\n'
    elif j == 7:
      line = line[:-1] + ' ' + str(max(dict_sample['x_L']))+'\n'
    elif j == 8:
      line = line[:-1] + ' ' + str(min(dict_sample['y_L']))+'\n'
    elif j == 9:
      line = line[:-1] + ' ' + str(max(dict_sample['y_L']))+'\n'
    elif j == 116:
      line = line[:-1] + "'"+str(dict_user['Mobility_eff'])+' '+str(dict_user['kappa_eta'])+" 1'\n"
    elif j == 138:
      line = line[:-1] + ' ' + str(dict_user['Energy_barrier'])+"'\n"
    elif j == 217:
      line = line[:-1] + ' ' + str(dict_user['dt_PF']*dict_user['n_t_PF']) +'\n'
    elif j == 221:
      line = line[:-1] + ' ' + str(dict_user['dt_PF']) +'\n'
    file_to_write.write(line)

  file_to_write.close()

# -----------------------------------------------------------------------------#

def sort_files_yade():
  '''
  Sort vtk files from yade simulation.
  '''
  # look for the next indice
  j = 0
  filepath = Path('vtk/2grains_'+str(j)+'.vtk')
  while filepath.exists():
      j = j + 1
      filepath = Path('vtk/2grains_'+str(j)+'.vtk')
  # rename
  os.rename('vtk/2grains-polyhedra-00000000.vtk','vtk/2grains_'+str(j)+'.vtk')
  os.rename('vtk/grain_1-polyhedra-00000000.vtk','vtk/grain1_'+str(j)+'.vtk')
  os.rename('vtk/grain_2-polyhedra-00000000.vtk','vtk/grain2_'+str(j)+'.vtk')
  os.rename('vtk/2grains-polyhedra-00000001.vtk','vtk/2grains_'+str(j+1)+'.vtk')
  os.rename('vtk/grain_1-polyhedra-00000001.vtk','vtk/grain1_'+str(j+1)+'.vtk')
  os.rename('vtk/grain_2-polyhedra-00000001.vtk','vtk/grain2_'+str(j+1)+'.vtk')

# -----------------------------------------------------------------------------#

def compare_volumes(dict_user, dict_sample):
    '''
    Compare the volume of the contact in Moose and in Yade.
    '''
    # Yade
    # already done in run_yade() in main.py

    # Moose
    # count
    counter = 0
    for i_y in range(len(dict_sample['y_L'])):
        for i_x in range(len(dict_sample['x_L'])):
            if dict_sample['eta_1_map'][i_y, i_x] > 0.5 and dict_sample['eta_2_map'][i_y, i_x] > 0.5:
                counter = counter + 1
    # adapt
    contact_volume_moose = counter * 1*(dict_sample['x_L'][1]-dict_sample['x_L'][0])*(dict_sample['y_L'][1]-dict_sample['y_L'][0])
    # tracker
    dict_user['L_contact_volume_moose'].append(contact_volume_moose)

    # plot
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    ax1.plot(dict_user['L_contact_volume_yade'], label='Yade')
    ax1.plot(dict_user['L_contact_volume_moose'], label='Moose')
    ax1.legend()
    ax1.set_title(r'Contact volume ($m^3$)',fontsize = 30)
    fig.savefig('plot/contact_volumes.png')
    plt.close(fig)

    # convert in nb of node
    L_nb_mesh_contact = []
    for i in range(len(dict_user['L_contact_volume_moose'])):
        L_nb_mesh_contact.append(dict_user['L_contact_volume_moose'][i]/(1*(dict_sample['x_L'][1]-dict_sample['x_L'][0])*(dict_sample['y_L'][1]-dict_sample['y_L'][0])))
    # plot
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    ax1.plot(L_nb_mesh_contact)
    ax1.set_title(r'Number of node in contact volume (-)',fontsize = 30)
    fig.savefig('plot/contact_nb_node.png')
    plt.close(fig)

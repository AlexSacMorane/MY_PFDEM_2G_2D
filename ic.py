# -*- encoding=utf-8 -*-

import numpy as np
import math

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def create_2_spheres(dict_user, dict_sample):
    '''
    Create initial conditions with 2 spheres.
    Mesh and phase field maps are generated
    '''
    # Create initial mesh
    print("Creating initial mesh")

    x_L = np.linspace(dict_user['x_min'], dict_user['x_max'], dict_user['n_mesh_x'])
    y_L = np.linspace(dict_user['y_min'], dict_user['y_max'], dict_user['n_mesh_y'])

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Create initial phase map

    print("Creating initial phase field maps")

    pos_1 = [0,-dict_user['radius']]
    pos_2 = [0, dict_user['radius']]

    eta_1_map = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    eta_2_map = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))

    # iteration on x
    for i_x in range(len(x_L)):
        x = x_L[i_x]
        # iteration on y
        for i_y in range(len(y_L)):
            y = y_L[i_y]

            # distance to g1 and g2
            d_node_to_g1 = np.linalg.norm(np.array([x,y])-np.array(pos_1))
            d_node_to_g2 = np.linalg.norm(np.array([x,y])-np.array(pos_2))

            # eta 1
            if d_node_to_g1 <= dict_user['radius']-dict_user['w_int']/2 :
                eta_1_map[-1-i_y, i_x] = 1
            elif dict_user['radius']-dict_user['w_int']/2 < d_node_to_g1 and d_node_to_g1 < dict_user['radius']+dict_user['w_int']/2:
                eta_1_map[-1-i_y, i_x] = 0.5*(1+math.cos(math.pi*(d_node_to_g1-dict_user['radius']+dict_user['w_int']/2)/dict_user['w_int']))
            elif dict_user['radius']+dict_user['w_int']/2 <= d_node_to_g1 :
                eta_1_map[-1-i_y, i_x] = 0

            # eta 2
            if d_node_to_g2 <= dict_user['radius']-dict_user['w_int']/2 :
                eta_2_map[-1-i_y, i_x] = 1
            elif dict_user['radius']-dict_user['w_int']/2 < d_node_to_g2 and d_node_to_g2 < dict_user['radius']+dict_user['w_int']/2:
                eta_2_map[-1-i_y, i_x] = 0.5*(1+math.cos(math.pi*(d_node_to_g2-dict_user['radius']+dict_user['w_int']/2)/dict_user['w_int']))
            elif dict_user['radius']+dict_user['w_int']/2 <= d_node_to_g2 :
                eta_2_map[-1-i_y, i_x] = 0

    # save dict
    dict_sample['pos_1'] = pos_1
    dict_sample['pos_2'] = pos_2
    dict_sample['eta_1_map'] = eta_1_map
    dict_sample['eta_2_map'] = eta_2_map
    dict_sample['x_L'] = x_L
    dict_sample['y_L'] = y_L

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def create_solute(dict_user, dict_sample):
    '''
    Create the map of the solute distribution.
    '''
    c_map = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    for i_x in range(len(dict_sample['x_L'])):
        for i_y in range(len(dict_sample['y_L'])):
            c_map[-1-i_y, i_x] = dict_user['C_eq'] # system at the equilibrium initialy
    # save in dict
    dict_sample['c_map'] = c_map

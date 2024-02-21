# -*- encoding=utf-8 -*-

import numpy as np
import math, skfmm

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
 
#------------------------------------------------------------------------------------------------------------------------------------------ #

def create_2_proxy_hexagons(dict_user, dict_sample):
    '''
    Create initial conditions with 2 hexagons like.
    Mesh and phase field maps are generated
    '''
    # Create initial mesh
    print("Creating initial mesh")

    x_L = np.linspace(dict_user['x_min'], dict_user['x_max'], dict_user['n_mesh_x'])
    y_L = np.linspace(dict_user['y_min'], dict_user['y_max'], dict_user['n_mesh_y'])

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Create initial phase map

    print("Creating initial phase field maps")

    pos_1 = [0,-dict_user['radius']*math.sqrt(3)/2]
    pos_2 = [0, dict_user['radius']*math.sqrt(3)/2]

    eta_1_map = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    eta_2_map = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))

    # iteration on x
    for i_x in range(len(x_L)):
        x = x_L[i_x]
        # iteration on y
        for i_y in range(len(y_L)):
            y = y_L[i_y]

            # g1 

            # distance
            d_node_to_g1 = np.linalg.norm(np.array([x,y])-np.array(pos_1))
            # orientation
            u_1 = np.array([x,y])-np.array(pos_1)
            if u_1[1] > 0 :
                phi_1 = math.acos(u_1[0]/d_node_to_g1)
            else :
                phi_1 = 2*math.pi - math.acos(u_1[0]/d_node_to_g1)
            # adapt orientation
            while not (0 <= phi_1 and phi_1 < math.pi/3):
                phi_1 = phi_1 - math.pi/3
            # compute the radius in the specific orientation
            if phi_1 < math.pi/6:
                r_phi_1 = dict_user['radius']*(1+(math.sqrt(3)-2)/2*phi_1/(math.pi/6)) 
            else :
                r_phi_1 = dict_user['radius']*(math.sqrt(3)/2+(2-math.sqrt(3))/2*(phi_1-math.pi/6)/(math.pi/6))
            # build eta1 in-out map
            if d_node_to_g1 < r_phi_1 : 
                eta_1_map[-1-i_y, i_x] = 0.5
            else :
                eta_1_map[-1-i_y, i_x] = -0.5

            # g2

            # distance    
            d_node_to_g2 = np.linalg.norm(np.array([x,y])-np.array(pos_2))
            # orientation
            u_2 = np.array([x,y])-np.array(pos_2)
            if u_2[1] > 0 :
                phi_2 = math.acos(u_2[0]/d_node_to_g2)
            else :
                phi_2 = 2*math.pi - math.acos(u_2[0]/d_node_to_g2)
            # adapt orientation
            while not (0 <= phi_2 and phi_2 < math.pi/3):
                phi_2 = phi_2 - math.pi/3
            # compute the radius in the specific orientation
            if phi_2 < math.pi/6:
                r_phi_2 = dict_user['radius']*(1+(math.sqrt(3)-2)/2*phi_2/(math.pi/6)) 
            else :
                r_phi_2 = dict_user['radius']*(math.sqrt(3)/2+(2-math.sqrt(3))/2*(phi_2-math.pi/6)/(math.pi/6))
            # build eta2 in-out map
            if d_node_to_g2 < r_phi_2 : 
                eta_2_map[-1-i_y, i_x] = 0.5
            else :
                eta_2_map[-1-i_y, i_x] = -0.5

    # compute sdf maps       
    sdf_1 = skfmm.distance(eta_1_map, dx = np.array([x_L[1]-x_L[0],y_L[1]-y_L[0]]))
    sdf_2 = skfmm.distance(eta_2_map, dx = np.array([x_L[1]-x_L[0],y_L[1]-y_L[0]]))

    # build phase maps
    # iteration on x
    for i_x in range(len(x_L)):
        # iteration on y
        for i_y in range(len(y_L)):

            # eta 1
            if sdf_1[i_y, i_x] > dict_user['w_int']/2 :
                eta_1_map[i_y, i_x] = 1
            elif sdf_1[i_y, i_x] < -dict_user['w_int']/2 :
                eta_1_map[i_y, i_x] = 0
            else :
                eta_1_map[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sdf_1[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

            # eta 2
            if sdf_2[i_y, i_x] > dict_user['w_int']/2 :
                eta_2_map[i_y, i_x] = 1
            elif sdf_2[i_y, i_x] < -dict_user['w_int']/2 :
                eta_2_map[i_y, i_x] = 0
            else :
                eta_2_map[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sdf_2[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

    # save dict
    dict_sample['pos_1'] = pos_1
    dict_sample['pos_2'] = pos_2
    dict_sample['eta_1_map'] = eta_1_map
    dict_sample['eta_2_map'] = eta_2_map
    dict_sample['x_L'] = x_L
    dict_sample['y_L'] = y_L
    
#------------------------------------------------------------------------------------------------------------------------------------------ #

def create_2_hexagons(dict_user, dict_sample):
    '''
    Create initial conditions with 2 hexagons.
    Mesh and phase field maps are generated
    '''
    # Create initial mesh
    print("Creating initial mesh")

    x_L = np.linspace(dict_user['x_min'], dict_user['x_max'], dict_user['n_mesh_x'])
    y_L = np.linspace(dict_user['y_min'], dict_user['y_max'], dict_user['n_mesh_y'])

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Create initial phase map

    print("Creating initial phase field maps")

    pos_1 = [0,-dict_user['radius']*math.sqrt(3)/2]
    pos_2 = [0, dict_user['radius']*math.sqrt(3)/2]

    eta_1_map = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))
    eta_2_map = np.zeros((dict_user['n_mesh_y'], dict_user['n_mesh_x']))

    # iteration on x
    for i_x in range(len(x_L)):
        x = x_L[i_x]
        # iteration on y
        for i_y in range(len(y_L)):
            y = y_L[i_y]

            # g1 

            # distance
            d_node_to_g1 = np.linalg.norm(np.array([x,y])-np.array(pos_1))
            # orientation
            u_1 = np.array([x,y])-np.array(pos_1)
            if u_1[1] > 0 :
                phi_1 = math.acos(u_1[0]/d_node_to_g1)
            else :
                phi_1 = 2*math.pi - math.acos(u_1[0]/d_node_to_g1)
            # adapt orientation
            while not (0 <= phi_1 and phi_1 < math.pi/3):
                phi_1 = phi_1 - math.pi/3
            # compute the radius in the specific orientation
            AC = dict_user['radius']
            C = math.pi/3 
            # build equation
            A_eq = AC**2*math.cos(phi_1)**2 - AC**2*math.cos(C)**2
            B_eq = 2*AC**3*math.cos(phi_1)*math.cos(C)**2 - 2*AC**3*math.cos(phi_1)
            C_eq = AC**4 - AC**4*math.cos(C)**2
            # solve equation
            Discri_eq = B_eq**2 - 4*A_eq*C_eq
            r_phi_1 = (-B_eq - math.sqrt(Discri_eq))/(2*A_eq)
            
            # build eta1 in-out map
            if d_node_to_g1 < r_phi_1 : 
                eta_1_map[-1-i_y, i_x] = 0.5
            else :
                eta_1_map[-1-i_y, i_x] = -0.5

            # g2

            # distance    
            d_node_to_g2 = np.linalg.norm(np.array([x,y])-np.array(pos_2))
            # orientation
            u_2 = np.array([x,y])-np.array(pos_2)
            if u_2[1] > 0 :
                phi_2 = math.acos(u_2[0]/d_node_to_g2)
            else :
                phi_2 = 2*math.pi - math.acos(u_2[0]/d_node_to_g2)
            # adapt orientation
            while not (0 <= phi_2 and phi_2 < math.pi/3):
                phi_2 = phi_2 - math.pi/3
            # compute the radius in the specific orientation
            AC = dict_user['radius']
            C = math.pi/3 
            # build equation
            A_eq = AC**2*math.cos(phi_2)**2 - AC**2*math.cos(C)**2
            B_eq = 2*AC**3*math.cos(phi_2)*math.cos(C)**2 - 2*AC**3*math.cos(phi_2)
            C_eq = AC**4 - AC**4*math.cos(C)**2
            # solve equation
            Discri_eq = B_eq**2 - 4*A_eq*C_eq
            r_phi_2 = (-B_eq - math.sqrt(Discri_eq))/(2*A_eq)
            # build eta2 in-out map
            if d_node_to_g2 < r_phi_2 : 
                eta_2_map[-1-i_y, i_x] = 0.5
            else :
                eta_2_map[-1-i_y, i_x] = -0.5

    # compute sdf maps       
    sdf_1 = skfmm.distance(eta_1_map, dx = np.array([x_L[1]-x_L[0],y_L[1]-y_L[0]]))
    sdf_2 = skfmm.distance(eta_2_map, dx = np.array([x_L[1]-x_L[0],y_L[1]-y_L[0]]))

    # build phase maps
    # iteration on x
    for i_x in range(len(x_L)):
        # iteration on y
        for i_y in range(len(y_L)):

            # eta 1
            if sdf_1[i_y, i_x] > dict_user['w_int']/2 :
                eta_1_map[i_y, i_x] = 1
            elif sdf_1[i_y, i_x] < -dict_user['w_int']/2 :
                eta_1_map[i_y, i_x] = 0
            else :
                eta_1_map[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sdf_1[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

            # eta 2
            if sdf_2[i_y, i_x] > dict_user['w_int']/2 :
                eta_2_map[i_y, i_x] = 1
            elif sdf_2[i_y, i_x] < -dict_user['w_int']/2 :
                eta_2_map[i_y, i_x] = 0
            else :
                eta_2_map[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sdf_2[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

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

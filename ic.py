# -*- encoding=utf-8 -*-

# -----------------------------------------------------------------------------#

def create_materials():
    '''
    Create materials.
    '''
    O.materials.append(FrictMat(young=-1, poisson=-1, frictionAngle=atan(0.5), density=2000, label='frictMat'))

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def create_2_spheres():
    '''
    Create initial conditions with 2 spheres.
    Potential blocks, mesh and phase field maps are generated
    '''
    # Create initial mesh
    print("Creating initial mesh")

    global x_L, y_L, z_L

    x_L = np.linspace(x_min, x_max, n_mesh_x)
    y_L = np.linspace(y_min, y_max, n_mesh_y)

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Create initial phase map

    print("Creating initial phase field maps")

    global pos_1, pos_2
    global eta_1_map, eta_2_map

    pos_1 = [0,-radius]
    pos_2 = [0, radius]

    eta_1_map = np.zeros((n_mesh_y, n_mesh_x))
    eta_2_map = np.zeros((n_mesh_y, n_mesh_x))

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
            if d_node_to_g1 <= radius-w/2 :
                eta_1_map[-1-i_y, i_x] = 1
            elif radius-w/2 < d_node_to_g1 and d_node_to_g1 < radius+w/2:
                eta_1_map[-1-i_y, i_x] = 0.5*(1+math.cos(math.pi*(d_node_to_g1-radius+w/2)/w))
            elif radius+w/2 <= d_node_to_g1 :
                eta_1_map[-1-i_y, i_x] = 0

            # eta 2
            if d_node_to_g2 <= radius-w/2 :
                eta_2_map[-1-i_y, i_x] = 1
            elif radius-w/2 < d_node_to_g2 and d_node_to_g2 < radius+w/2:
                eta_2_map[-1-i_y, i_x] = 0.5*(1+math.cos(math.pi*(d_node_to_g2-radius+w/2)/w))
            elif  radius+w/2 <= d_node_to_g2 :
                eta_2_map[-1-i_y, i_x] = 0

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # Interpolate

    global L_n_plane_1, L_d_plane_1, L_n_plane_2, L_d_plane_2
    L_n_plane_1, L_d_plane_1 = interpolate_planes(eta_1_map, pos_1)
    L_n_plane_2, L_d_plane_2 = interpolate_planes(eta_2_map, pos_2)

    # ------------------------------------------------------------------------------------------------------------------------------------------ #
    # create potential blocks

    create_potential_block()

# ------------------------------------------------------------------------------------------------------------------------------------------ #

def create_solute():
    '''
    Create the map of the solute distribution.
    '''
    global c_map
    c_map = np.zeros((n_mesh_y, n_mesh_x))
    for i_x in range(len(x_L)):
        for i_y in range(len(y_L)):
            c_map[-1-i_y, i_x] = C_eq # system at the equilibrium initialy

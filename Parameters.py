#------------------------------------------------------------------------------------------------------------------------------------------ #
# Librairies
#------------------------------------------------------------------------------------------------------------------------------------------#

import numpy as np

#------------------------------------------------------------------------------------------------------------------------------------------ #
# Parameters
#------------------------------------------------------------------------------------------------------------------------------------------#

def get_parameters():
    '''
    Define the parameters used in the simulation.
    '''

    #---------------------------------------------------------------------#
    # PFDEM

    n_DEMPF_ite = 50 # number of PFDEM iterations
    n_proc = 5 # number of processors used
    j_total = 0 # index global of results
    save_simulation = False # indicate if the simulation is saved
    n_max_vtk_files = 10 # maximum number of vtk files (can be None to save all files)

    # Select Figures to plot
    # Available:
    # contact_pressure, contact_distrib_m_ed, contact_point_ed, contact_volume, contact_nb_node
    # contact_detection, contact_h_s_v, contact_dem
    # m_ed, dt_PF, IC, processor, sphericities, force_applied, maps, n_grain_kc_map, dim_dom
    # shape_evolution, n_vertices, sum_etai_c, mean_etai_c, mass_loss, performances
    # disp_strain_andrade, sample_height, y_contactPoint
    L_figures = ['maps', 'shape_evolution', 'dt_PF', 'IC', 'dim_dom', 'contact_volume',\
                 'contact_h_s_v', 'contact_distrib_m_ed', 'contact_point_ed', 'disp_strain_andrade']

    # Figure (plot all or current)
    # The maps configuration
    print_all_map_config = False # else only the current one is printed

    #---------------------------------------------------------------------#
    # DEM (Yade)

    # steady state detection
    n_ite_max = 5000 # maximum number of iteration during a DEM step
    n_steady_state_detection = 100 # number of iterations considered in the window
    # the difference between max and min < tolerance * force_applied
    steady_state_detection = 0.02
    # + the force applied must be contained in this window

    # sollicitation
    force_applied = 3e6 # N
    control_force = True # Boolean to determine if the force is controled with the contact volume

    # DEM material parameters
    # Young modulus
    E = 1e8 # Pa
    # Poisson ratio
    Poisson = 0.3

    # Figure (plot all or current)
    # The evolution of the overlap durig DEM steps
    print_all_contact_dem = False # else only the current one is printed
    # The evolution of shape (compared to the initial state)
    print_all_shape_evolution = False # else only the current one is printed

    #---------------------------------------------------------------------#
    # Grain description

    # shape of the grain
    # Sphere, Hex or Proxy_Hex
    Shape = 'Sphere'

    # the radius of grains
    radius = 1 # m
    # discretization of the grain
    n_phi = 60

    #---------------------------------------------------------------------#
    # Phase-Field (Moose)

    # determine if remeshing is available
    remesh = True
    # mesh
    check_database = True
    # remesh
    if remesh:
        size_x_mesh = 0.02
        size_y_mesh = size_x_mesh
        m_size_mesh = (size_x_mesh+size_y_mesh)/2
        margin_mesh_domain = 10*size_x_mesh
    # constant mesh
    else :
        x_min = -1.3*radius
        x_max =  1.3*radius
        y_min = -2.3*radius
        y_max =  2.3*radius
        n_mesh_x = 100
        n_mesh_y = 200
        m_size_mesh = ((x_max-x_min)/(n_mesh_x-1)+(y_max-y_min)/(n_mesh_y-1))/2

    # PF material parameters
    # the energy barrier
    Energy_barrier = 1
    # number of mesh in the interface
    n_int = 6
    # the interface thickness
    w = m_size_mesh*n_int
    # the gradient coefficient
    kappa_eta = Energy_barrier*w*w/9.86
    # the mobility
    Mobility_eff = 1

    # kinetics of dissolution and precipitation
    # it affects the tilting coefficient in Ed
    k_diss = 0.1 # mol.m-2.s-1
    k_prec = 0.1

    # molar concentration at the equilibrium
    C_eq = 1 # number of C_ref, mol m-3

    # diffusion of the solute
    D_solute = 10 # m2 s-1
    n_struct_element = int(round(radius*0.10/m_size_mesh,0))
    struct_element = np.array(np.ones((n_struct_element,n_struct_element)), dtype=bool) # for dilation

    # Aitken method
    # the time stepping and duration of one PF simualtion
    # level 0
    dt_PF_0 = 0.2 # time step
    # level 1
    dt_PF_1 = dt_PF_0/2
    m_ed_contact_1 = 0.1
    # level 2
    dt_PF_2 = dt_PF_1/2
    m_ed_contact_2 = 0.2
    # n_t_PF*dt_PF gives the total time duration
    n_t_PF = 10 # number of iterations

    # Contact box detection
    eta_contact_box_detection = 0.25 # value of the phase field searched to determine the contact box

    # Figure (plot all or current)
    # The detection of the contact by a box
    print_all_contact_detection = False # else only the current one is printed

    #---------------------------------------------------------------------#
    # trackers

    L_displacement = []
    L_sum_eta_1 = []
    L_sum_eta_2 = []
    L_sum_c = []
    L_sum_mass = []
    L_m_eta_1 = []
    L_m_eta_2 = []
    L_m_c = []
    L_m_mass = []
    L_distance_extrema = []
    L_equivalent_area = []
    L_contact_overlap = []
    L_contact_area = []
    L_contact_volume_yade = []
    L_contact_volume_moose = []
    L_contact_volume_box = []
    L_t_pf_to_dem_1 = []
    L_t_pf_to_dem_2 = []
    L_t_dem = []
    L_t_dem_to_pf = []
    L_t_pf = []
    L_P_applied = []
    L_n_v_1 = []
    L_n_v_2 = []
    L_n_v_1_target = []
    L_n_v_2_target = []
    L_m_ed = []
    L_m_ed_contact = []
    L_m_ed_large_contact = []
    L_m_ed_plus_contact = []
    L_m_ed_minus_contact = []
    L_m_ed_plus_large_contact = []
    L_m_ed_minus_large_contact = []
    L_ed_contact_point = []
    L_vertices_1_init = None
    L_force_applied = []
    L_dt_PF = []
    L_AreaSphericity = []
    L_DiameterSphericity = []
    L_CircleRatioSphericity = []
    L_PerimeterSphericity = []
    L_WidthToLengthRatioSpericity = []
    L_grain_kc_map = []
    L_sample_height = []
    L_y_contactPoint = []
    L_loss_move_pf_eta1 = []
    L_loss_move_pf_eta2 = []
    L_loss_move_pf_c = []
    L_loss_move_pf_m = []
    L_loss_kc_eta1 = []
    L_loss_kc_eta2 = []
    L_loss_kc_c = []
    L_loss_kc_m = []
    L_loss_pf_eta1 = []
    L_loss_pf_eta2 = []
    L_loss_pf_c = []
    L_loss_pf_m = []
    if remesh:
        L_x_min_dom = []
        L_x_max_dom = []
        L_y_min_dom = []
        L_y_max_dom = []
        L_delta_x_max = []
        L_delta_y_max = []

    #---------------------------------------------------------------------#
    # dictionnary

    dict_user = {
    'n_DEMPF_ite': n_DEMPF_ite,
    'n_proc': n_proc,
    'j_total': j_total,
    'save_simulation': save_simulation,
    'n_max_vtk_files': n_max_vtk_files,
    'L_figures': L_figures,
    'print_all_map_config': print_all_map_config,
    'n_ite_max': n_ite_max,
    'n_steady_state_detection': n_steady_state_detection,
    'steady_state_detection': steady_state_detection,
    'force_applied': force_applied,
    'force_applied_target': force_applied,
    'control_force': control_force,
    'E': E,
    'Poisson': Poisson,
    'print_all_contact_dem': print_all_contact_dem,
    'print_all_shape_evolution': print_all_shape_evolution,
    'Shape': Shape,
    'radius': radius,
    'remesh': remesh,
    'n_phi': n_phi,
    'Mobility_eff': Mobility_eff,
    'kappa_eta': kappa_eta,
    'n_int': n_int,
    'w_int': w,
    'Energy_barrier': Energy_barrier,
    'dt_PF_0': dt_PF_0,
    'dt_PF_1': dt_PF_1,
    'Aitken_1': m_ed_contact_1,
    'dt_PF_2': dt_PF_2,
    'Aitken_2': m_ed_contact_2,
    'n_t_PF': n_t_PF,
    'k_diss': k_diss,
    'k_prec': k_prec,
    'C_eq': C_eq,
    'D_solute': D_solute,
    'struct_element': struct_element,
    'eta_contact_box_detection': eta_contact_box_detection,
    'print_all_contact_detection': print_all_contact_detection,
    'L_displacement': L_displacement,
    'L_sum_eta_1': L_sum_eta_1,
    'L_sum_eta_2': L_sum_eta_2,
    'L_sum_c': L_sum_c,
    'L_sum_mass': L_sum_mass,
    'L_m_eta_1': L_m_eta_1,
    'L_m_eta_2': L_m_eta_2,
    'L_m_c': L_m_c,
    'L_m_mass': L_m_mass,
    'L_distance_extrema': L_distance_extrema,
    'L_equivalent_area': L_equivalent_area,
    'L_contact_overlap': L_contact_overlap,
    'L_contact_area': L_contact_area,
    'L_contact_volume_yade': L_contact_volume_yade,
    'L_contact_volume_moose': L_contact_volume_moose,
    'L_contact_volume_box': L_contact_volume_box,
    'L_t_pf_to_dem_1': L_t_pf_to_dem_1,
    'L_t_pf_to_dem_2': L_t_pf_to_dem_2,
    'L_t_dem': L_t_dem,
    'L_t_dem_to_pf': L_t_dem_to_pf,
    'L_t_pf': L_t_pf,
    'L_P_applied': L_P_applied,
    'L_n_v_1': L_n_v_1,
    'L_n_v_2': L_n_v_2,
    'L_n_v_1_target': L_n_v_1_target,
    'L_n_v_2_target': L_n_v_2_target,
    'L_m_ed': L_m_ed,
    'L_m_ed_contact': L_m_ed_contact,
    'L_m_ed_large_contact': L_m_ed_large_contact,
    'L_m_ed_plus_contact': L_m_ed_plus_contact,
    'L_m_ed_minus_contact': L_m_ed_minus_contact,
    'L_m_ed_plus_large_contact': L_m_ed_plus_large_contact,
    'L_m_ed_minus_large_contact': L_m_ed_minus_large_contact,
    'L_ed_contact_point': L_ed_contact_point,
    'L_vertices_1_init': L_vertices_1_init,
    'L_force_applied': L_force_applied,
    'L_dt_PF': L_dt_PF,
    'L_AreaSphericity': L_AreaSphericity,
    'L_DiameterSphericity': L_DiameterSphericity,
    'L_CircleRatioSphericity': L_CircleRatioSphericity,
    'L_PerimeterSphericity': L_PerimeterSphericity,
    'L_WidthToLengthRatioSpericity': L_WidthToLengthRatioSpericity,
    'L_grain_kc_map': L_grain_kc_map,
    'L_sample_height': L_sample_height,
    'L_y_contactPoint': L_y_contactPoint,
    'L_loss_move_pf_eta1': L_loss_move_pf_eta1,
    'L_loss_move_pf_eta2': L_loss_move_pf_eta2,
    'L_loss_move_pf_c': L_loss_move_pf_c,
    'L_loss_move_pf_m': L_loss_move_pf_m,
    'L_loss_kc_eta1': L_loss_kc_eta1,
    'L_loss_kc_eta2': L_loss_kc_eta2,
    'L_loss_kc_c': L_loss_kc_c,
    'L_loss_kc_m': L_loss_kc_m,
    'L_loss_pf_eta1': L_loss_pf_eta1,
    'L_loss_pf_eta2': L_loss_pf_eta2,
    'L_loss_pf_c': L_loss_pf_c,
    'L_loss_pf_m': L_loss_pf_m,
    'check_database': check_database
    }

    # specific inputs
    if remesh:
        dict_user['size_x_mesh'] = size_x_mesh
        dict_user['size_y_mesh'] = size_y_mesh
        dict_user['margin_mesh_domain'] = margin_mesh_domain
        dict_user['L_x_min_dom'] = L_x_min_dom
        dict_user['L_x_max_dom'] = L_x_max_dom
        dict_user['L_y_min_dom'] = L_y_min_dom
        dict_user['L_y_max_dom'] = L_y_max_dom
        dict_user['L_delta_x_max'] = L_delta_x_max
        dict_user['L_delta_y_max'] = L_delta_y_max
    else :
        dict_user['x_min'] = x_min
        dict_user['x_max'] = x_max
        dict_user['y_min'] = y_min
        dict_user['y_max'] = y_max
        dict_user['n_mesh_x'] = n_mesh_x
        dict_user['n_mesh_y'] = n_mesh_y

    return dict_user

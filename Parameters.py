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

    n_DEMPF_ite = 100 # number of PFDEM iterations
    n_proc = 6 # number of processors used
    j_total = 0 # index global of results

    #---------------------------------------------------------------------#
    # DEM (Yade)

    # steady state detection
    n_ite_max = 2000 # maximum number of iteration during a DEM step
    n_steady_state_detection = 100 # number of iterations considered in the window
    # the difference between max and min < tolerance * force_applied
    steady_state_detection = 0.02
    # + the force applied must be contained in this window

    # sollicitation
    force_applied = 2e7 # N

    # DEM material parameters
    # Young modulus
    E = 1e9 # Pa
    # Poisson ratio
    Poisson = 0.3

    #---------------------------------------------------------------------#
    # Grain description

    # the radius of grains
    radius = 1 # m
    # discretization of the grain
    n_phi = 40

    #---------------------------------------------------------------------#
    # Phase-Field (Moose)

    # mesh
    x_min = -1.15*radius
    x_max =  1.15*radius
    n_mesh_x = 100
    y_min = -2.15*radius
    y_max =  2.15*radius
    n_mesh_y = 200

    # PF material parameters
    # the energy barrier
    Energy_barrier = 1
    # number of mesh in the interface
    n_int = 6
    # the interface thickness
    w = ((x_max-x_min)/(n_mesh_x-1)+(y_max-y_min)/(n_mesh_y-1))*n_int
    # the gradient coefficient
    kappa_eta = Energy_barrier*w*w/9.68
    # the mobility
    Mobility = 3/2.2*w
    Mobility_eff = 2.2/3*Mobility/w
    # the time stepping and duration of one PF simualtion
    dt_PF = 0.0005 # time step
    n_t_PF = 5 # number of iterations
    # n_t_PF*dt_PF gives the total time duration

    # kinetics of dissolution and precipitation
    # it affects the tilting coefficient in Ed
    k_diss = 0.5 # mol.m-2.s-1
    k_prec = k_diss

    # molar concentration at the equilibrium
    C_eq = 1 # number of C_ref, mol m-3

    # diffusion of the solute
    D_solute = 1000 # m2 s-1
    struct_element = np.array(np.ones((3,3)), dtype=bool) # for dilation

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

    #---------------------------------------------------------------------#
    # dictionnary

    dict_user = {
    'n_DEMPF_ite': n_DEMPF_ite,
    'n_proc': n_proc,
    'j_total': j_total,
    'n_ite_max': n_ite_max,
    'n_steady_state_detection': n_steady_state_detection,
    'steady_state_detection': steady_state_detection,
    'force_applied': force_applied,
    'E': E,
    'Poisson': Poisson,
    'radius': radius,
    'n_phi': n_phi,
    'x_min': x_min,
    'x_max': x_max,
    'n_mesh_x': n_mesh_x,
    'y_min': y_min,
    'y_max': y_max,
    'n_mesh_y': n_mesh_y,
    'Mobility': Mobility,
    'Mobility_eff': Mobility_eff,
    'kappa_eta': kappa_eta,
    'n_int': n_int,
    'w_int': w,
    'Energy_barrier': Energy_barrier,
    'dt_PF': dt_PF,
    'n_t_PF': n_t_PF,
    'k_diss': k_diss,
    'k_prec': k_prec,
    'C_eq': C_eq,
    'D_solute': D_solute,
    'struct_element': struct_element,
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
    'L_m_ed': L_m_ed
    }

    return dict_user

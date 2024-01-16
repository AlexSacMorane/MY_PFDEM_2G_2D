#------------------------------------------------------------------------------------------------------------------------------------------ #
# Parameters
#------------------------------------------------------------------------------------------------------------------------------------------#

def get_parameters():
    '''
    Define the parameters used in the simulation.
    '''

    #---------------------------------------------------------------------#
    # PFDEM

    n_DEMPF_ite = 20 # number of PFDEM iterations
    n_proc = 6 # number of processors used
    j_total = 0 # index global of results

    #---------------------------------------------------------------------#
    # DEM (Yade)

    # steady state detection
    n_ite_max = 10000 # maximum number of iteration during a DEM step
    n_steady_state_detection = 100 # number of iterations considered in the window
    # the difference between max and min < tolerance * force_applied
    steady_state_detection = 0.01
    # + the force applied must be contained in this window

    # sollicitation
    force_applied = 1e2

    # DEM material parameters
    # Normal stiffness
    Kn = 1e3  #Pa/m
    # ?? stiffness
    Ks = Kn * 2 / 3  #Pa/m

    #---------------------------------------------------------------------#
    # Grain description

    # the radius of grains
    radius = 1 # m
    r = 0.1 # inner particule
    n_phi = 10 # discretization of the grain

    #---------------------------------------------------------------------#
    # Phase-Field (Moose)

    # mesh
    x_min = -1.1*radius
    x_max =  1.1*radius
    n_mesh_x = 50
    y_min = -2.1*radius
    y_max =  2.1*radius
    n_mesh_y = 100

    # PF material parameters
    # the mobility
    Mobility = 1
    # the gradient coefficient
    kappa_eta = 0.02
    # the interface thickness
    w = ((x_max-x_min)/2/(n_mesh_x-1)+(y_max-y_min)/2/(n_mesh_y-1))*5
    # the enrgy barrier
    Energy_barrier = 20*kappa_eta/w/w

    # the time stepping and duration of one PF simualtion
    dt_PF = 2 # time step
    n_t_PF = 1 # number of iterations
    # n_t_PF*dt_PF gives the total time duration

    # kinetics of dissolution and precipitation
    # it affects the tilting coefficient in Ed
    k_diss = 1 # mol.m-2.s-1
    k_prec = k_diss

    # molar concentration at the equilibrium
    C_eq = 1 # number of C_ref, mol m-3

    # diffusion of the solute
    D_solute = 100 # m2 s-1

    #---------------------------------------------------------------------#
    # trackers

    L_displacement = []
    L_sum_eta_1 = []
    L_sum_eta_2 = []
    L_sum_c = []

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
    'Kn': Kn,
    'Ks': Ks,
    'radius': radius,
    'r_pot': r,
    'n_phi': n_phi,
    'x_min': x_min,
    'x_max': x_max,
    'n_mesh_x': n_mesh_x,
    'y_min': y_min,
    'y_max': y_max,
    'n_mesh_y': n_mesh_y,
    'Mobility': Mobility,
    'kappa_eta': kappa_eta,
    'w_int': w,
    'Energy_barrier': Energy_barrier,
    'dt_PF': dt_PF,
    'n_t_PF': n_t_PF,
    'k_diss': k_diss,
    'k_prec': k_prec,
    'C_eq': C_eq,
    'D_solute': D_solute,
    'L_displacement': L_displacement,
    'L_sum_eta_1': L_sum_eta_1,
    'L_sum_eta_2': L_sum_eta_2,
    'L_sum_c': L_sum_c
    }

    return dict_user

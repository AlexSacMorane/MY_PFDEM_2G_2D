# -*- encoding=utf-8 -*-

# -----------------------------------------------------------------------------#

def move_phasefield():
    '''
    Move phase field maps by interpolation.
    '''
    # g2 (as g1 is fixed)
    displacement_g2 = O.bodies[1].state.pos - O.bodies[1].state.refPos
    L_displacement.append(displacement_g2[2])

    # loading old variables
    global eta_1_map, eta_2_map, c_map

    print('Updating phase field maps')

    # updating phase map
    eta_1_map_new = np.zeros((n_mesh_y, n_mesh_x))
    eta_2_map_new = np.zeros((n_mesh_y, n_mesh_x))
    c_map_new = c_map # not updated

    # iteration on y
    i_y_old = 0
    for i_y in range(len(y_L)):
        y = y_L[i_y]

        # eta 1, fixed
        eta_1_map_new = eta_1_map

        # eta 2
        if displacement_g2[1] < 0:
            if y-displacement_g2[1] <= y_L[-1]:
                # look for window
                while not (y_L[i_y_old] <= y-displacement_g2[1] and y-displacement_g2[1] < y_L[i_y_old+1]):
                    i_y_old = i_y_old + 1
                # interpolate
                eta_2_map_new[-1-i_y, :] = (eta_2_map[-1-(i_y_old+1), :] - eta_2_map[-1-i_y_old, :])/(y_L[i_y_old+1] - y_L[i_y_old])*\
                                           (y-displacement_g2[1] - y_L[i_y_old]) + eta_2_map[-1-i_y_old, :]
        elif displacement_g2[1] > 0:
            if y_L[0] <= y-displacement_g2[1]:
                # look for window
                while not (y_L[i_y_old] <= y-displacement_g2[1] and y-displacement_g2[1] < y_L[i_y_old+1]):
                    i_y_old = i_y_old + 1
                # interpolate
                eta_2_map_new[-1-i_y, :] = (eta_2_map[-1-(i_y_old+1), :] - eta_2_map[-1-i_y_old, :])/(y_L[i_y_old+1] - y_L[i_y_old])*\
                                           (y-displacement_g2[1] - y_L[i_y_old]) + eta_2_map[-1-i_y_old, :]
        else :
            eta_2_map_new = eta_2_map

        # solute
        # TO DO : push solute with grain displacement

    # update variables
    eta_1_map = eta_1_map_new
    eta_2_map = eta_2_map_new
    c_map = c_map_new

    # write txts for phase field
    write_eta_txt() # phase field
    write_c_txt() # solute

# -----------------------------------------------------------------------------#

def compute_ed():
    '''
    Compute sink/source term.
    '''
    global ed_map
    ed_map = np.zeros((n_mesh_y,n_mesh_x))

    contact_area = O.interactions[0,1].phys.contactArea
    R_gas = 82.06e5 # cm3 Pa K-1 mol-1
    Temp = 25+278   # K
    V_m = 27.1      # cm3 mol-1
    conv = 1        #  = 1/C_ref, normalize the sink/source
    v_mesh = (x_L[1]-x_L[0])*(y_L[1]-y_L[0]) # surface of a mesh, m2

    for i_x in range(len(x_L)):
        for i_y in range(len(y_L)):
            # determine pressure
            if eta_1_map[i_y, i_x] > 0.5 and eta_2_map[i_y, i_x] > 0.5: # in the contact
                P =  np.linalg.norm(O.interactions[0,1].phys.normalForce)/contact_area # Pa
            else : # not in the contact
                P = 0 # Pa
            r = contact_area*math.exp(P*V_m/(R_gas*Temp))*(1-c_map[i_y, i_x]/(C_eq*math.exp(P*V_m/(R_gas*Temp))))
            if r >= 0 :
                r = k_diss*r # dissolution rate
            else :
                r = k_prec*r # precipitation rate
            # the rate is in mol s-1
            r = r/v_mesh # convert in mol m-3 s-1
            r = r*conv # convert in C_ref s-1
            # save in the map
            ed_map[i_y, i_x] = r

    # write ed
    write_ed_txt()

# -----------------------------------------------------------------------------#

def compute_ep_old():
    '''
    Compute mechanical term.
    '''
    ep_map = np.zeros((n_mesh, n_mesh, n_mesh))

    r_map = np.zeros((n_mesh, n_mesh, n_mesh))

    contact_area = O.interactions[0,1].phys.contactArea
    Ah = 0.5114 # T = 25°C
    Bh = 0.3288 # T = 25°C
    R_const = 82.057 # pas constant in cm3 atm K-1 mol-1
    Temp_C = 25 # °C
    Temp_K = Temp_C+273 # K
    # Solubility at 1 bar
    K_1T = 6.13 # mol kg water
    # Halite
    V_Halite = 27.1 *10**-6 # m3 mol-1
    # Cl
    a_o_Cl = 3.63 * 10**-10 # m
    bi_Cl = 0.017 # -
    z_Cl = -1 # -
    a_1_Cl =  44.65    # cal bar-1 mol-1
    a_2_Cl =  4.801e-2 # cal mol-1
    a_3_Cl =  4.325    # cal K-1 bar-1 mol-1
    a_4_Cl = -2.847e-4 # cal-1 K-1 mol-1
    w_Cl = 1.748e-5    # cal-1 mol-1
    b_1_Cl = -0.331 # -
    b_2_Cl =  20.16 # -
    b_3_Cl =  0     # -
    b_4_Cl =  1     # -
    # Na
    a_o_Na = 4.08 * 10**-10 # m
    bi_Na = 0.082 # -
    z_Na =  1 # -
    a_1_Na =  22.8     # cal bar-1 mol-1
    a_2_Na = -4.38e-2  # cal mol-1
    a_3_Na = -4.10     # cal K-1 bar-1 mol-1
    a_4_Na = -0.586e-4 # cal-1 K-1 mol-1
    w_Na = 0.09e-5     # cal-1 mol-1
    b_1_Na = 0.3       # -
    b_2_Na = 52        # -
    b_3_Na = -3.33e-3  # -
    b_4_Na = 0.566     # -
    # dielectric constant
    U_1 =  3.4279e2  # -
    U_2 = -5.0866e-3 # °C-1
    U_3 =  9.4690e-7 # °C-2
    U_4 = -2.0525    # -
    U_5 =  3.159e3   # °C
    U_6 = -1.8289e2  # °C
    U_7 = -8.0325e3  # bar
    U_8 =  4.2142e6  # bar °C
    U_9 =  2.1417    # bar / °C
    alpha    = (U_7 + U_8/Temp_C + U_9*Temp_C) # bar
    beta     = U_4 + U_5/(U_6+Temp_C)
    eps_1000 = U_1*math.exp(U_2*Temp_C+U_3*Temp_C**2)
    # compressibility of pure water
    kappa_0 = 4.52e-5 # bar-1
    # Debye-Hückel
    A_gamma = 0.51 #mol kg-1 water a t 25°C
    B_gamma = 1 # Debye-Hückel lenght parameter (1/A_o)

    for i_x in range(n_mesh):
        for i_y in range(n_mesh):
            for i_z in range(n_mesh):
                # determine pressure
                if eta_1_map[i_x, i_y, i_z] > 0.5 and eta_2_map[i_x, i_y, i_z] > 0.5: # in the contact
                    P =  O.interactions[0,1].phys.normalForce/contact_area*1e-5 # bar
                else : # not in the contact
                    P = 1 # bar
                # ionic strenght
                I = 0.5*c_map[i_x,i_y,i_z]*z_Cl**2+0.5*c_map[i_x,i_y,i_z]*z_Na**2 # c must be in mol/l
                # activity coefficient
                gamma_Na = 10**(-Ah*z_Na**2*I**0.5/(1+a_o_Na*Bh*I**0.5)+bi_Na*I)
                gamma_Cl = 10**(-Ah*z_Cl**2*I**0.5/(1+a_o_Cl*Bh*I**0.5)+bi_Cl*I)
                # ion activity product
                Q = gamma_Cl*gamma_Na*c_map[i_x,i_y,i_z]**2 # diffusion c_Na and c_Cl are the same
                # derivatives of dielectric constant
                d_eps_inv = beta/((alpha+P)*(beta*math.log((alpha+P)/(alpha+1e3))+eps_1000)**2)
                d_ln_eps  = beta/((alpha+P)*(beta*math.log((alpha+P)/(alpha+1e3))+eps_1000))
                # Debye-Hückel limiting slope
                Av = 1.534*A_gamma*R_const*Temp_K*(3*d_ln_eps-kappa_0)
                # molar volume
                Vm_Na_0 = 41.84*(0.1*a_1_Na+100*a_2_Na/(2600+P)+a_3_Na/(Temp_K-228)+10**4*a_4_Na/((2600+P)*(Temp_K-228))-w_Na*d_eps_inv)
                Vm_Na = Vm_Na_0 + 0.5*Av*z_Na**2*I**0.5/(1+a_o_Na*B_gamma*I**0.5) + (b_1_Na + b_2_Na/(Temp_K-228) + b_3_Na*(Temp_K-228))*I**b_4_Na
                Vm_Cl_0 = 41.84*(0.1*a_1_Cl+100*a_2_Cl/(2600+P)+a_3_Cl/(Temp_K-228)+10**4*a_4_Cl/((2600+P)*(Temp_K-228))-w_Cl*d_eps_inv)
                Vm_Cl = Vm_Cl_0 + 0.5*Av*z_Cl**2*I**0.5/(1+a_o_Cl*B_gamma*I**0.5) + (b_1_Cl + b_2_Cl/(Temp_K-228) + b_3_Cl*(Temp_K-228))*I**b_4_Cl
                # volume change of reaction
                Delta_Vr = Vm_Na + Vm_Cl - V_Halite
                # solubility
                Keq = K_1T*10**(-Delta_Vr*(P-1)/R_const*Temp_K)
                # sink/source term
                r_map[i_x, i_y, i_z] = - contact_area*k_diss*(1-Q/Keq)

    write_ep_txt(ep_map)

#-------------------------------------------------------------------------------

def write_eta_txt():
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
    for x in x_L:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    file_to_write_2.write(line)
    # y
    file_to_write_1.write('AXIS Y\n')
    file_to_write_2.write('AXIS Y\n')
    line = ''
    for y in y_L:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    file_to_write_2.write(line)
    # data
    file_to_write_1.write('DATA\n')
    file_to_write_2.write('DATA\n')
    for j in range(len(y_L)):
        for i in range(len(x_L)):
            # grain 1
            file_to_write_1.write(str(eta_1_map[-1-j,i])+'\n')
            # grain 2
            file_to_write_2.write(str(eta_2_map[-1-j,i])+'\n')
    # close
    file_to_write_1.close()
    file_to_write_2.close()

#-------------------------------------------------------------------------------

def write_c_txt():
    '''
    Write a .txt file needed for MOOSE simulation.

    This .txt represents the solute map.
    '''
    file_to_write = open('data/c.txt','w')
    # x
    file_to_write.write('AXIS X\n')
    line = ''
    for x in x_L:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # y
    file_to_write.write('AXIS Y\n')
    line = ''
    for y in y_L:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # data
    file_to_write.write('DATA\n')
    for j in range(len(y_L)):
        for i in range(len(x_L)):
            file_to_write.write(str(c_map[-1-j,i])+'\n')
    # close
    file_to_write.close()

#-------------------------------------------------------------------------------

def write_ed_txt():
    '''
    Write a .txt file needed for MOOSE simulation.

    This .txt represents the map of the external mechanical energy.
    '''
    file_to_write = open('data/ed.txt','w')
    # x
    file_to_write.write('AXIS X\n')
    line = ''
    for x in x_L:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # y
    file_to_write.write('AXIS Y\n')
    line = ''
    for y in y_L:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write.write(line)
    # data
    file_to_write.write('DATA\n')
    for j in range(len(y_L)):
        for i in range(len(x_L)):
            file_to_write.write(str(ed_map[-1-j,i])+'\n')
    # close
    file_to_write.close()

#-------------------------------------------------------------------------------

def compute_kc():
    '''
    Compute the diffusion coefficient of the solute.
    Then Write a .txt file needed for MOOSE simulation.

    This .txt file represent the phase field maps.
    '''
    # compute
    kc_map = np.array(np.zeros((n_mesh_y, n_mesh_x)))
    # iterate on x and y
    for i_y in range(len(y_L)):
        for i_x in range(len(x_L)):
            if eta_1_map[i_y, i_x] < 0.5 and eta_2_map[i_y, i_x] < 0.5: # out of the grain
                kc_map[i_y, i_x] = D_solute
            elif eta_1_map[i_y, i_x] > 0.5 and eta_2_map[i_y, i_x] > 0.5: # in the contact
                kc_map[i_y, i_x] = D_solute
            else :
                kc_map[i_y, i_x] = 0

    # write
    file_to_write_1 = open('data/kc.txt','w')
    # x
    file_to_write_1.write('AXIS X\n')
    line = ''
    for x in x_L:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    # y
    file_to_write_1.write('AXIS Y\n')
    line = ''
    for y in y_L:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write_1.write(line)
    # data
    file_to_write_1.write('DATA\n')
    for j in range(len(y_L)):
        for i in range(len(x_L)):
            # grain 1
            file_to_write_1.write(str(kc_map[-1-j,i])+'\n')
    # close
    file_to_write_1.close()


#-------------------------------------------------------------------------------

def write_i():
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
      line = line[:-1] + ' ' + str(len(x_L)-1)+'\n'
    elif j == 5:
      line = line[:-1] + ' ' + str(len(y_L)-1)+'\n'
    elif j == 6:
      line = line[:-1] + ' ' + str(min(x_L))+'\n'
    elif j == 7:
      line = line[:-1] + ' ' + str(max(x_L))+'\n'
    elif j == 8:
      line = line[:-1] + ' ' + str(min(y_L))+'\n'
    elif j == 9:
      line = line[:-1] + ' ' + str(max(y_L))+'\n'
    elif j == 116:
      line = line[:-1] + "'"+str(Mobility)+' '+str(kappa_eta)+" 1'\n"
    elif j == 138:
      line = line[:-1] + ' ' + str(Energy_barrier)+"'\n"
    elif j == 217:
      line = line[:-1] + ' ' + str(dt_PF*n_t_PF) +'\n'
    elif j == 221:
      line = line[:-1] + ' ' + str(dt_PF) +'\n'
    file_to_write.write(line)

  file_to_write.close()

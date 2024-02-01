# -*- encoding=utf-8 -*-

import pickle, math, os, shutil
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import vtk
from vtk.util.numpy_support import vtk_to_numpy

# -----------------------------------------------------------------------------#

def sort_files(dict_user, dict_sample):
     '''
     Sort files generated by MOOSE to different directories
     '''
     os.rename('pf_out.e','vtk/pf_out.e')
     os.rename('pf.i','input/pf.i')
     j = 0
     j_str = index_to_str(j)
     j_total_str = index_to_str(dict_user['j_total'])
     filepath = Path('pf_other_'+j_str+'.pvtu')
     while filepath.exists():
         for i_proc in range(dict_user['n_proc']):
            os.rename('pf_other_'+j_str+'_'+str(i_proc)+'.vtu','vtk/pf_other_'+j_str+'_'+str(i_proc)+'.vtu')
            shutil.copyfile('vtk/pf_other_'+j_str+'_'+str(i_proc)+'.vtu','vtk/pf_'+j_total_str+'_'+str(i_proc)+'.vtu')
         os.rename('pf_other_'+j_str+'.pvtu','vtk/pf_other_'+j_str+'.pvtu')
         # write .pvtu to save all vtk
         file = open('vtk/pf_'+j_total_str+'.pvtu','w')
         file.write('''<?xml version="1.0"?>
         <VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">
         \t<PUnstructuredGrid GhostLevel="1">
         \t\t<PPointData>
         \t\t\t<PDataArray type="Float64" Name="ed"/>
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
             line = line + '''\t\t<Piece Source="pf_'''+j_total_str+'''_'''+str(i_proc)+'''.vtu"/>\n'''
         file.write(line)
         file.write('''\t</PUnstructuredGrid>
         </VTKFile>''')
         file.close()
         j = j + 1
         j_str = index_to_str(j)
         filepath = Path('pf_other_'+j_str+'.pvtu')
         dict_user['j_total'] = dict_user['j_total'] + 1
         j_total_str = index_to_str(dict_user['j_total'])
     return index_to_str(j-1)

# -----------------------------------------------------------------------------#

def index_to_str(j):
    '''
    Convert a integer into a string with the format XXX.
    '''
    if j < 10:
        return '00'+str(j)
    elif 10<=j and j<100:
        return '0'+str(j)
    else :
        return str(j)

# -----------------------------------------------------------------------------#

def read_vtk(dict_user, dict_sample, j_str):
    '''
    Read the last vtk files to obtain data from MOOSE.

    Do not work calling yade.
    '''
    eta_1_map_old = dict_sample['eta_1_map'].copy()
    eta_2_map_old = dict_sample['eta_2_map'].copy()
    c_map_old = dict_sample['c_map'].copy()
    L_XYZ = []
    L_eta1 = []
    L_eta2 = []
    L_c = []
    L_limits = []

    # iterate on the proccessors used
    for i_proc in range(dict_user['n_proc']):

        # name of the file to load
        namefile = 'vtk/pf_other_'+j_str+'_'+str(i_proc)+'.vtu'

        # load a vtk file as input
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(namefile)
        reader.Update()

        # Grab a scalar from the vtk file
        nodes_vtk_array = reader.GetOutput().GetPoints().GetData()
        eta1_vtk_array = reader.GetOutput().GetPointData().GetArray("eta1")
        eta2_vtk_array = reader.GetOutput().GetPointData().GetArray("eta2")
        c_vtk_array = reader.GetOutput().GetPointData().GetArray("c")

        #Get the coordinates of the nodes and the scalar values
        nodes_array = vtk_to_numpy(nodes_vtk_array)
        eta1_array = vtk_to_numpy(eta1_vtk_array)
        eta2_array = vtk_to_numpy(eta2_vtk_array)
        c_array = vtk_to_numpy(c_vtk_array)

        # look for limits
        x_min = None
        x_max = None
        y_min = None
        y_max = None

        # Must detect common zones between processors
        for i_XYZ in range(len(nodes_array)) :
            XYZ = nodes_array[i_XYZ]
            # Do not consider twice a point
            if list(XYZ) not in L_XYZ :
                L_XYZ.append(list(XYZ))
                L_eta1.append(eta1_array[i_XYZ])
                L_eta2.append(eta2_array[i_XYZ])
                L_c.append(c_array[i_XYZ])
                # set first point
                if x_min == None :
                    x_min = list(XYZ)[0]
                    x_max = list(XYZ)[0]
                    y_min = list(XYZ)[1]
                    y_max = list(XYZ)[1]
                # look for limits of the processor
                else :
                    if list(XYZ)[0] < x_min:
                        x_min = list(XYZ)[0]
                    if list(XYZ)[0] > x_max:
                        x_max = list(XYZ)[0]
                    if list(XYZ)[1] < y_min:
                        y_min = list(XYZ)[1]
                    if list(XYZ)[1] > y_max:
                        y_max = list(XYZ)[1]

        # Here the algorithm can be help as the mapping is known

        # save limits
        L_limits.append([x_min,x_max,y_min,y_max])

    # plot processors distribution
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    # parameters
    title_fontsize = 20
    for i_proc in range(len(L_limits)):
        limits = L_limits[i_proc]
        ax1.plot([limits[0],limits[1],limits[1],limits[0],limits[0]],[limits[2],limits[2],limits[3],limits[3],limits[2]], label='proc '+str(i_proc))
    ax1.legend()
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_title('Processor i has the priority on i+1',fontsize = title_fontsize)
    fig.suptitle('Processors ditribution',fontsize = 1.2*title_fontsize)
    fig.savefig('plot/processors_distribution.png')
    plt.close(fig)

    # rebuild map from lists
    for i_XYZ in range(len(L_XYZ)):
        # find nearest node
        L_search = list(abs(np.array(dict_sample['x_L']-L_XYZ[i_XYZ][0])))
        i_x = L_search.index(min(L_search))
        L_search = list(abs(np.array(dict_sample['y_L']-L_XYZ[i_XYZ][1])))
        i_y = L_search.index(min(L_search))
        # rewrite map
        dict_sample['eta_1_map'][-1-i_y, i_x] = L_eta1[i_XYZ]
        dict_sample['eta_2_map'][-1-i_y, i_x] = L_eta2[i_XYZ]
        dict_sample['c_map'][-1-i_y, i_x] = L_c[i_XYZ]

# -----------------------------------------------------------------------------#

def read_vtk_own(j_str):
    '''
    Read the last vtk files to obtain data from MOOSE.
    '''
    global eta_1_map, eta_2_map, c_map
    eta_1_map_old = eta_1_map.copy()
    eta_2_map_old = eta_2_map.copy()
    c_map_old = c_map.copy()

    id_L = None
    data_jump_len = len('          ')

    for i_proc in range(n_proc):
        L_Work = [[], #X
                  [], #Y
                  [], #eta1
                  [], #eta2
                  []] #c

        file = open('vtk/pf_other_'+j_str+'_'+str(i_proc)+'.vtu','r')
        lines = file.readlines()
        file.close()
        #iterations on line
        for line in lines:

            if '"eta1"' in line:
                id_L = 2
            if '"eta2"' in line:
                id_L = 3
            if '"c"' in line:
                id_L = 4
            if '"Points"' in line:
                id_L = 0
            if ('</DataArray>' in line or  '<InformationKey' in line) and id_L != None:
                id_L = None
            if line[0:data_jump_len] == '          ' and (id_L == 2 or id_L == 3): #Read etai
                line = line[data_jump_len:]
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        L_Work[id_L].append(float(line[c_start:c_end]))
                        c_start = c_i+1
                L_Work[id_L].append(float(line[c_start:]))
            if line[0:data_jump_len] == '          ' and id_L == 4: #Read c
                line = line[data_jump_len:]
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        L_Work[id_L].append(float(line[c_start:c_end]))
                        c_start = c_i+1
                L_Work[id_L].append(float(line[c_start:]))

            if line[0:data_jump_len] == '          ' and id_L == 0: #Read [X, Y]
                line = line[data_jump_len:]
                XYZ_temp = []
                c_start = 0
                for c_i in range(0,len(line)):
                    if line[c_i]==' ':
                        c_end = c_i
                        XYZ_temp.append(float(line[c_start:c_end]))
                        if len(XYZ_temp)==3:
                            L_Work[0].append(XYZ_temp[0])
                            L_Work[1].append(XYZ_temp[1])
                            XYZ_temp = []
                        c_start = c_i+1
                XYZ_temp.append(float(line[c_start:]))
                L_Work[0].append(XYZ_temp[0])
                L_Work[1].append(XYZ_temp[1])
        #Adaptating data
        for i in range(len(L_Work[0])):
            #Interpolation method
            L_dx = []
            for x_i in x_L :
                L_dx.append(abs(x_i - L_Work[0][i]))
            L_dy = []
            for y_i in y_L :
                L_dy.append(abs(y_i - L_Work[1][i]))
            # update phase field maps
            eta_1_map[-1-list(L_dy).index(min(L_dy)), list(L_dx).index(min(L_dx))] = L_Work[2][i]
            eta_2_map[-1-list(L_dy).index(min(L_dy)), list(L_dx).index(min(L_dx))] = L_Work[3][i]
            # update solute map
            c_map[-1-list(L_dy).index(min(L_dy)), list(L_dx).index(min(L_dx))] = L_Work[4][i]

    # interpolate plane where etai = 0.5
    global L_n_plane_1, L_d_plane_1, L_n_plane_2, L_d_plane_2
    L_n_plane_1, L_d_plane_1 = interpolate_planes(eta_1_map, pos_1)
    L_n_plane_2, L_d_plane_2 = interpolate_planes(eta_2_map, pos_2)

# -----------------------------------------------------------------------------#

def compute_plane(dict_user, dict_sample):
    '''
    From a phase map, compute planes for the potential block
    '''
    # compute planes
    L_n_plane_1, L_d_plane_1 = interpolate_planes(dict_sample['eta_1_map'], dict_sample['pos_1'], dict_user, dict_sample)
    L_n_plane_2, L_d_plane_2 = interpolate_planes(dict_sample['eta_2_map'], dict_sample['pos_2'], dict_user, dict_sample)

    # save data
    dict_save = {
    'L_n_plane_1': L_n_plane_1,
    'L_d_plane_1': L_d_plane_1,
    'L_n_plane_2': L_n_plane_2,
    'L_d_plane_2': L_d_plane_2
    }
    with open('data/planes.data', 'wb') as handle:
        pickle.dump(dict_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

# -----------------------------------------------------------------------------#

def interpolate_planes(eta_i_map, center, dict_user, dict_sample):
    '''
    Interpolate plane for potential block.
    '''
    map_phi = []
    L_phi = []
    for i_phi in range(dict_user['n_phi']):
        phi = 2*math.pi*i_phi/dict_user['n_phi']
        L_phi.append(phi)
        map_phi.append([])
    L_phi.append(2*math.pi)
    for i_x in range(len(dict_sample['x_L'])-1):
        for i_y in range(len(dict_sample['y_L'])-1):
            L_in = [] # list the nodes inside the grain
            if eta_i_map[-1-i_y    , i_x] > 0.5 :
                L_in.append(0)
            if eta_i_map[-1-(i_y+1), i_x] > 0.5 :
                L_in.append(1)
            if eta_i_map[-1-(i_y+1), i_x+1] > 0.5 :
                L_in.append(2)
            if eta_i_map[-1-i_y    , i_x+1] > 0.5 :
                L_in.append(3)
            if L_in != [] and L_in != [0,1,2,3]:
                center_mesh = (np.array([dict_sample['x_L'][i_x], dict_sample['y_L'][i_y]])+np.array([dict_sample['x_L'][i_x+1], dict_sample['y_L'][i_y+1]]))/2
                u = (center_mesh-np.array(center))/np.linalg.norm(center_mesh-np.array(center))
                # compute phi
                if u[1]>=0:
                    phi = math.acos(u[0])
                else :
                    phi = 2*math.pi-math.acos(u[0])
                # iterate on the lines of the mesh to find the plane intersection
                L_p = []
                if (0 in L_in and 1 not in L_in) or (0 not in L_in and 1 in L_in):# line 01
                    x_p = dict_sample['x_L'][i_x]
                    y_p = (0.5-eta_i_map[-1-i_y, i_x])/(eta_i_map[-1-(i_y+1), i_x]-eta_i_map[-1-i_y, i_x])*(dict_sample['y_L'][i_y+1]-dict_sample['y_L'][i_y])+dict_sample['y_L'][i_y]
                    L_p.append(np.array([x_p, y_p])-np.array(center))
                if (1 in L_in and 2 not in L_in) or (1 not in L_in and 2 in L_in):# line 12
                    x_p = (0.5-eta_i_map[-1-(i_y+1), i_x])/(eta_i_map[-1-(i_y+1), i_x+1]-eta_i_map[-1-(i_y+1), i_x])*(dict_sample['x_L'][i_x+1]-dict_sample['x_L'][i_x])+dict_sample['x_L'][i_x]
                    y_p = dict_sample['y_L'][i_y+1]
                    L_p.append(np.array([x_p, y_p])-np.array(center))
                if (2 in L_in and 3 not in L_in) or (2 not in L_in and 3 in L_in):# line 23
                    x_p = dict_sample['x_L'][i_x+1]
                    y_p = (0.5-eta_i_map[-1-i_y, i_x+1])/(eta_i_map[-1-(i_y+1), i_x+1]-eta_i_map[-1-i_y, i_x+1])*(dict_sample['y_L'][i_y+1]-dict_sample['y_L'][i_y])+dict_sample['y_L'][i_y]
                    L_p.append(np.array([x_p, y_p])-np.array(center))
                if (3 in L_in and 0 not in L_in) or (3 not in L_in and 0 in L_in):# line 30
                    x_p = (0.5-eta_i_map[-1-i_y, i_x])/(eta_i_map[-1-i_y, i_x+1]-eta_i_map[-1-i_y, i_x])*(dict_sample['x_L'][i_x+1]-dict_sample['x_L'][i_x])+dict_sample['x_L'][i_x]
                    y_p = dict_sample['y_L'][i_y]
                    L_p.append(np.array([x_p, y_p])-np.array(center))
                # compute the mean point
                p_mean = np.array([0,0])
                for p in L_p :
                    p_mean = p_mean + p
                p_mean = p_mean/len(L_p)
                # look phi in L_phi
                i_phi = 0
                while not (L_phi[i_phi] <= phi and phi < L_phi[i_phi+1]) :
                    i_phi = i_phi + 1
                # save p_mean in the map
                map_phi[i_phi].append(p_mean)

    L_n_plane = []
    L_d_plane = []
    # interpolate plane (Least squares method)
    for i_phi in range(len(map_phi)):
        # compute phi
        phi = 2*math.pi*i_phi/dict_user['n_phi']
        # verify if a y=ax+b or a x=ay+b is better
        if math.sin(phi) > 0.5 or math.sin(phi) < -0.5:
            # compute plane equation
            result_ls = least_square_y_fx(map_phi[i_phi])
            slope = result_ls[0]
            cst = result_ls[1]
            # compute normal
            n_plane = np.array([-slope, 1])/np.linalg.norm(np.array([-slope, 1]))
            if math.sin(phi)*n_plane[1] < 0 :
                n_plane = -n_plane
            # compute distance
            d_plane = abs(cst)/math.sqrt(slope**2+(-1)**2)

        else :
            # compute plane equation
            result_ls = least_square_x_fy(map_phi[i_phi])
            slope = result_ls[0]
            cst = result_ls[1]
            # compute normal
            n_plane = np.array([-1, slope])/np.linalg.norm(np.array([-1, slope]))
            if math.sin(phi)*n_plane[1] < 0 :
                n_plane = -n_plane
            # compute distance
            d_plane = abs(cst)/math.sqrt((-1)**2+slope**2)

        # save
        L_n_plane.append(n_plane)
        L_d_plane.append(d_plane)

    return L_n_plane, L_d_plane

# -----------------------------------------------------------------------------#

def least_square_y_fx(L_p):
    '''
    y = result[0]*x + result[1]
    '''
    Matrix = np.array(np.zeros((2,2)))
    Vector = np.array(np.zeros(2))
    # compute terms
    for p in L_p:
        Matrix[0,0] = Matrix[0,0] + p[0]*p[0]
        Matrix[0,1] = Matrix[0,1] + p[0]
        Matrix[1,0] = Matrix[1,0] + p[0]
        Matrix[1,1] = Matrix[1,1] + 1
        Vector = Vector + np.array([p[0]*p[1], p[1]])
    # solve
    result = np.linalg.solve(Matrix, Vector)
    return result

# -----------------------------------------------------------------------------#

def least_square_x_fy(L_p):
    '''
    x = result[0]*y + result[1]
    '''
    Matrix = np.array(np.zeros((2,2)))
    Vector = np.array(np.zeros(2))
    # compute terms
    for p in L_p:
        Matrix[0,0] = Matrix[0,0] + p[1]*p[1]
        Matrix[0,1] = Matrix[0,1] + p[1]
        Matrix[1,0] = Matrix[1,0] + p[1]
        Matrix[1,1] = Matrix[1,1] + 1
        Vector = Vector + np.array([p[1]*p[0], p[0]])
    # solve
    result = np.linalg.solve(Matrix, Vector)
    return result

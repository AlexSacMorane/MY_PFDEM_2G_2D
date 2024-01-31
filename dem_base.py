# -*- encoding=utf-8 -*-

from yade import pack, utils, plot, export
import polyhedra_utils
import pickle
import numpy as np

# Own functions
execfile('dem.py')

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Load data

# from main
with open('data/main_to_dem.data', 'rb') as handle:
    dict_save = pickle.load(handle)
E = dict_save['E']
Poisson = dict_save['Poisson']
force_applied = dict_save['force_applied']
pos_1 = dict_save['pos_1']
pos_2 = dict_save['pos_2']
n_ite_max = dict_save['n_ite_max']
steady_state_detection = dict_save['steady_state_detection']
n_steady_state_detection = dict_save['n_steady_state_detection']
i_DEMPF_ite = dict_save['i_DEMPF_ite']

# from plane interpolation
with open('data/planes.data', 'rb') as handle:
    dict_save = pickle.load(handle)
L_vertices_1 = dict_save['L_vertices_1']
L_vertices_2 = dict_save['L_vertices_2']

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Plan simulation

# vtk exporter
vtkExporter = export.VTKExporter('vtk/2grains')
vtkExporter_1 = export.VTKExporter('vtk/grain_1')
vtkExporter_2 = export.VTKExporter('vtk/grain_2')

# Plot
create_plots() # from dem.py
# materials
create_materials() # from dem.py
# create grains
create_polyhedral() # from dem.py
# Engines
create_engines() # from dem.py
# time step
compute_dt() # from dem.py

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# DEM

O.run()
O.wait()

# ------------------------------------------------------------------------------------------------------------------------------------------ #
# Output

displacement_g2 = np.array(O.bodies[1].state.pos - O.bodies[1].state.refPos)
contact_volume = float(O.interactions[0,1].geom.penetrationVolume)
contact_area = float(O.interactions[0,1].geom.equivalentCrossSection) # Volume**(2/3)
contact_overlap = float(O.interactions[0,1].geom.equivalentPenetrationDepth) # Volume/Surface
normalForce = np.array(O.interactions[0,1].phys.normalForce)
pos_1 = [float(O.bodies[0].state.pos[0]), float(O.bodies[0].state.pos[1])]
pos_2 = [float(O.bodies[1].state.pos[0]), float(O.bodies[1].state.pos[1])]
overlap = 0
for v_i in O.bodies[0].shape.v:
    for v_j in O.bodies[1].shape.v:
        overlap_try = (O.bodies[0].state.pos[1]+v_i[1]) - (O.bodies[1].state.pos[1]+v_j[1])
        if overlap_try > overlap:
            overlap = overlap_try
surface = contact_volume/overlap
n_v_1 = len(O.bodies[0].shape.v)
n_v_2 = len(O.bodies[1].shape.v)

# Save data
dict_save = {
'displacement_g2': displacement_g2,
'contact_area': contact_area,
'contact_volume': contact_volume,
'contact_overlap': contact_overlap,
'distance_extrema': overlap,
'equivalent_area': surface,
'normalForce': normalForce,
'pos_1': pos_1,
'pos_2': pos_2,
'n_v_1': n_v_1,
'n_v_2': n_v_2,
'n_v_1_target': len(L_vertices_1),
'n_v_2_target': len(L_vertices_2)
}
with open('data/dem_to_main.data', 'wb') as handle:
    pickle.dump(dict_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

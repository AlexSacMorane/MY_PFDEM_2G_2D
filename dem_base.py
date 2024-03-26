# -*- encoding=utf-8 -*-

from yade import pack, utils, plot, export
import polyhedra_utils
import pickle
import numpy as np

# Own functions
# -----------------------------------------------------------------------------#

def create_materials():
    '''
    Create materials.
    '''
    O.materials.append(PolyhedraMat(young=E, poisson=Poisson, frictionAngle=atan(0.5), density=2000, label='frictMat'))

# -----------------------------------------------------------------------------#

def create_polyhedral():
    '''
    Recreate polydra from data extrapolated with phase field output.
    '''
    print("Creating polyhedra")

    # g1, bottom grain, fixed
    O.bodies.append(
        polyhedra_utils.polyhedra(
            O.materials[-1],
            v = L_vertices_1,
            fixed = True,
            color = [0,0,1]
        )
    )
    O.bodies[-1].state.refPos = O.bodies[-1].state.pos

    # g2, top grain
    O.bodies.append(
        polyhedra_utils.polyhedra(
            O.materials[-1],
            v = L_vertices_2,
            fixed = True,
            color = [1,0,0]
        )
    )
    O.bodies[-1].state.blockedDOFs = 'xzXYZ'
    O.bodies[-1].state.refPos = O.bodies[-1].state.pos

    # initial export
    vtkExporter.exportPolyhedra(what=dict(color='b.shape.color'))
    vtkExporter_1.exportPolyhedra(ids=[0])
    vtkExporter_2.exportPolyhedra(ids=[1])

# -----------------------------------------------------------------------------#

def applied_force():
    '''
    Apply force on the top grain
    '''
    O.forces.addF(1, (0, -force_applied, 0))

# -----------------------------------------------------------------------------#

def add_data():
    '''
    Add data to plot :
        - iteration
        - overlap between the two particles (>0 if overlap)
    '''
    if O.interactions.has(0,1) : # check if interaction exists
        if O.interactions[0,1].isReal: # check if interaction is real
            overlap = O.interactions[0,1].geom.equivalentPenetrationDepth
            volume = O.interactions[0,1].geom.penetrationVolume
            normal_force = O.interactions[0,1].phys.normalForce[1]
        else :
            overlap = 0
            volume = 0
            normal_force = 0
    else :
        overlap = 0
        volume = 0
        normal_force = 0
    plot.addData(iteration=O.iter, overlap=overlap, volume=volume, normal_force=normal_force, force_applied=force_applied)

# -----------------------------------------------------------------------------#

def check():
    '''
    Try to detect a wteady-state of the overlap between the two particles.
    A maximum number of iteration is used.
    '''
    if O.iter < max(n_ite_max*0.01, n_steady_state_detection):
        return
    window = plot.data['normal_force'][-n_steady_state_detection:]
    if O.iter > n_ite_max or \
       ((max(window)-min(window))<steady_state_detection*force_applied and
        max(window)>force_applied and min(window)<force_applied):
        vtkExporter.exportPolyhedra(what=dict(color='b.shape.color')) # final export
        vtkExporter_1.exportPolyhedra(ids=[0])
        vtkExporter_2.exportPolyhedra(ids=[1])
        if print_contact_dem:
            if print_all_contact_dem:
                plot.plot(noShow=True).savefig('plot/contact_dem/'+str(i_DEMPF_ite)+'.png')
            else:
                plot.plot(noShow=True).savefig('plot/contact_dem.png')
        O.pause() # stop DEM simulation

# -----------------------------------------------------------------------------#

def compute_dt():
    '''
    Compute the time step used in the DEM step.
    '''
    O.dt = 0.1 * polyhedra_utils.PWaveTimeStep()

# -----------------------------------------------------------------------------#

def create_engines():
    '''
    Create engines.

    Ip2_PolyhedraMat_PolyhedraMat_PolyhedraPhys 
    Normal: 1/kn = 1/Y1 + 1/Y2 
    Shear: 1/ks = 1/Y1v1 + 1/Y2v2
    Y is the Young modulus
    v is the Poisson ratio

    Law2_PolyhedraGeom_PolyhedraPhys_Volumetric 
    F = k N
    Force is proportionnal to the volume
    '''
    O.engines = [
            ForceResetter(),
            InsertionSortCollider([Bo1_Polyhedra_Aabb()], verletDist=0.00),
            InteractionLoop(
                    [Ig2_Polyhedra_Polyhedra_PolyhedraGeom()],
                    [Ip2_PolyhedraMat_PolyhedraMat_PolyhedraPhys()],
                    [Law2_PolyhedraGeom_PolyhedraPhys_Volumetric()]
            ),
    		PyRunner(command='applied_force()',iterPeriod=1),
            NewtonIntegrator(damping=0.5, exactAsphericalRot=True, gravity=[0, 0, 0], label='newton'),
    		PyRunner(command='add_data()',iterPeriod=1),
            PyRunner(command='check()',iterPeriod=1, label='checker')
    ]

# -----------------------------------------------------------------------------#

def create_plots():
    '''
    Create plots during the DEM step.
    '''
    plot.plots = {'iteration': ('overlap', None, 'volume'),'iteration ':('normal_force','force_applied')}

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
print_all_contact_dem = dict_save['print_all_contact_dem']
print_contact_dem = dict_save['print_contact_dem']
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
'n_v_2_target': len(L_vertices_2),
'contact_point': np.array(O.interactions[0,1].geom.contactPoint)
}
with open('data/dem_to_main.data', 'wb') as handle:
    pickle.dump(dict_save, handle, protocol=pickle.HIGHEST_PROTOCOL)

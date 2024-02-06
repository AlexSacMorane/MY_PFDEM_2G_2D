# -*- encoding=utf-8 -*-

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

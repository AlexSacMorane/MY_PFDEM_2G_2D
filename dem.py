# -*- encoding=utf-8 -*-

# -----------------------------------------------------------------------------#

def create_materials():
    '''
    Create materials.
    '''
    O.materials.append(FrictMat(young=-1, poisson=-1, frictionAngle=atan(0.5), density=2000, label='frictMat'))

# -----------------------------------------------------------------------------#

def create_potential_block():
    '''
    Recreate potential blocks from data extrapolated with phase field output.
    '''
    print("Creating potential block")

    # g1, bottom grain, fixed
    a_L = []
    b_L = []
    c_L = []
    d_L = []
    for i in range(len(L_n_plane_1)):
        a_L.append(float(L_n_plane_1[i][0]))
        b_L.append(float(L_n_plane_1[i][1]))
        c_L.append(0)
        d_L.append(float(abs(L_d_plane_1[i]))-r)
    # z+
    a_L.append(0)
    b_L.append(0)
    c_L.append(1)
    d_L.append(2*r-r)
    # z-
    a_L.append(0)
    b_L.append(0)
    c_L.append(-1)
    d_L.append(2*r-r)
    g1 = Body()
    g1.aspherical = True
    g1.shape = PotentialBlock(
            a=a_L,
            b=b_L,
            c=c_L,
            d=d_L,
            r=r,
            R=0.0,
            AabbMinMax=True,
            isBoundary=True,
            color = [0,0,1]
    )
    utils._commonBodySetup(g1, g1.shape.volume, g1.shape.inertia, material='frictMat', pos=[pos_1[0], pos_1[1], 0], fixed=True)
    g1.state.ori = g1.shape.orientation
    O.bodies.append(g1)

    # g2, top grain
    a_L = []
    b_L = []
    c_L = []
    d_L = []
    for i in range(len(L_n_plane_2)):
        a_L.append(float(L_n_plane_2[i][0]))
        b_L.append(float(L_n_plane_2[i][1]))
        c_L.append(0)
        d_L.append(float(abs(L_d_plane_2[i]))-r)
    # z+
    a_L.append(0)
    b_L.append(0)
    c_L.append(1)
    d_L.append(2*r-r)
    # z-
    a_L.append(0)
    b_L.append(0)
    c_L.append(-1)
    d_L.append(2*r-r)
    g2 = Body()
    g2.aspherical = True
    g2.shape = PotentialBlock(
            a=a_L,
            b=b_L,
            c=c_L,
            d=d_L,
            r=r,
            R=0.0,
            AabbMinMax=True,
            isBoundary=False,
            color = [1,0,0]
    )
    utils._commonBodySetup(g2, g2.shape.volume, g2.shape.inertia, material='frictMat', pos=[pos_2[0], pos_2[1], 0], fixed=True, blockedDOFs='xzXYZ')
    g2.state.ori = g2.shape.orientation
    O.bodies.append(g2)

    # initial export
    vtkExporter.exportPotentialBlocks(what=dict(color='b.shape.color'))
    vtkExporter_1.exportPotentialBlocks(ids=[0])
    vtkExporter_2.exportPotentialBlocks(ids=[1])

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
            overlap = O.interactions[0,1].geom.penetrationDepth
            volume = O.interactions[0,1].phys.contactArea*overlap
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
        vtkExporter.exportPotentialBlocks(what=dict(color='b.shape.color')) # final export
        vtkExporter_1.exportPotentialBlocks(ids=[0])
        vtkExporter_2.exportPotentialBlocks(ids=[1])
        plot.plot(noShow=True).savefig('plot/contact/contact_'+str(i_DEMPF_ite)+'.png')
        O.pause() # stop DEM simulation

# -----------------------------------------------------------------------------#

def compute_dt():
    '''
    Compute the time step used in the DEM step.
    '''
    O.dt = 0.1 * sqrt(min(O.bodies[0].state.mass, O.bodies[1].state.mass) / Kn)

# -----------------------------------------------------------------------------#

def create_engines():
    '''
    Create engines.
    '''
    O.engines = [
            ForceResetter(),
            InsertionSortCollider([PotentialBlock2AABB()], verletDist=0.00),
            InteractionLoop(
                    [Ig2_PB_PB_ScGeom(calContactArea=True)],
                    [Ip2_FrictMat_FrictMat_KnKsPBPhys(kn_i=Kn, ks_i=Ks, Knormal=Kn, Kshear=Ks, viscousDamping=0.2)],
                    [Law2_SCG_KnKsPBPhys_KnKsPBLaw(label='law', initialOverlapDistance = 1e-6, allowViscousAttraction=False)]
            ),
    		PyRunner(command='applied_force()',iterPeriod=1),
            NewtonIntegrator(damping=0.0, exactAsphericalRot=True, gravity=[0, 0, 0], label='newton'),
    		PyRunner(command='add_data()',iterPeriod=1),
            PyRunner(command='check()',iterPeriod=1, label='checker')
    ]

# -----------------------------------------------------------------------------#

def create_plots():
    '''
    Create plots during the DEM step.
    '''
    plot.plots = {'iteration': ('overlap', None, 'volume'),'iteration ':('normal_force','force_applied')}

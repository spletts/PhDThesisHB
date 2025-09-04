"""Run very similar simulation setup as in simCRpropa https://github.com/me-manu/simCRpropa, 
but with the purpose of understanding the trajectory of a charged particle in magnetic field with different max step sizes.
"""

from crpropa import *
import os
import numpy as np


# Track all particles
THINNING = 0
# Origin for source, observer sphere, magnetic field grid
ORIGIN =  Vector3d(0, 0, 0)
# If True, track particle trajectory (rather than collect events on sphere)
TRAJECTORY = True
# Particle ID. Inject electron (charged particle to study impact of step size and magnetic field on electrons/positrons)
PID = int(11)
# Magnitude of magnetic field at grid points in Gauss. Note the grid is interpolated to solve for mangetic field at arbitrary positions
BAMP = 1e-13
# Inject 1 particle (no need to track more than 1 particle trajectory)
NINJ = 1
# Redshift. Note this impacts the photon background fields and magnetic field etc due to comoving coordinates and/or evolution 
REDSHIFTS = [0.1]  # [0.1, 0.4, 0.7]
# Max step sizes in Mpc. For 1 TeV electron, the Larmor radius is about 1e-2
MAX_STEP_MPC = [50, 1, 1e-6]  # [50, 10, 5, 1, 0.1, 1e-6]
# Units: eV. 110 MeV is just above breaking condition of 1e8 eV. The highest energy *photons* injected ar 1e13 eV (~10 TeV), may split energy in half via pair production
E_INJ = 1.1e8
# Limits for magnetic field plot, in Mpc
B_AXLIMS = 2e-5
# Max trajectory to trace in Mpc
# This can be less than mean free path (not interested in interactions, just charged particle spiraling in B field), and 
# also small enough to resolve the Larmor radius (long trajectory distance can swamp out the tight spiral in a plot), or make a cut in `plot_trajectory.py` to resolve spiral (~2e-5 Mpc)
MAX_TRAJ_MPC = [1]  # [5, 1, 0.1]  # [0.001, 1]
# Tolerance (units: kpc?)
TOL = 1e-9


def create_bfield(grid_len_meters, bamp, cell_len_mpc=50, gridOrigin=ORIGIN, axlims_mpc=None, outdir=None, random_direction=False):
    """Setup magnetic field with CRPropa modules

    Parameters
    ----------
    grid_len_meters : float
        Size of each side of base grid (not individual grid cells). Units of meters.
    gridOrigin : CRPropa Vector3d
        Origin for magnetic field grid
    random_direction : bool
        If True use random numbers to define direction of magnetic field.
        If False, used fixed direction at each grid point as in `initFixedDirectionField`

    Remainder defined in `simcrpropa_main` docstring.
    """
    # Determine size of each side of grid cell in order to have `ngrid` cells between the source and observer
    ngrid = int(np.ceil(grid_len_meters / (cell_len_mpc * Mpc)))  # dsrc has units of meters
    print(f"grid spacing: {cell_len_mpc} Mpc ({cell_len_mpc*Mpc} meters) {ngrid} cells")
    # B field grid
    print(ngrid, cell_len_mpc*Mpc)
    vgrid = Grid3f(gridOrigin,
                    ngrid,
                    cell_len_mpc*Mpc)  # meters
    if random_direction:
        initRandomField(vgrid, Bamplitude=bamp*gauss)
    else:
        initFixedDirectionField(vgrid, Bamplitude=bamp*gauss)
    bfield0 = MagneticFieldGrid(vgrid)
    # field, extent, origin, reflective
    bfield = PeriodicMagneticField(bfield0, Vector3d(ngrid*cell_len_mpc*Mpc), gridOrigin, True)

    if axlims_mpc is not None:
        plotBField(bfield, axlims_mpc, outdir)

    return bfield


def create_source(redshift, pid, inj_powlaw, delta_einj_list=None, src_loc=ORIGIN):
    """Setup CRPropa source modules.

    Parameters
    ----------
    src_loc : CRPropa Vector3d
        Put the source at this location (meters).

    Remainder defined in `simcrpropa_main` docstring.
    """
    source = Source()
    source.add(SourcePosition(src_loc))
    source.add(SourceRedshift(redshift))
    source.add(SourceParticleType(pid))

    source.add(SourceEmissionCone(Vector3d(1, 1, 0), np.deg2rad(0.1)))

    if inj_powlaw:
        # TODO update this function if using powerlaw
        source.add(SourcePowerLawSpectrum(4 * GeV, 6 * GeV, -1.5))
    else:
        for delta_einj in delta_einj_list:
            source.add(SourceEnergy(delta_einj*eV))

    print(source)

    return source


def get_outdir(pid, num_inj, redshift, powlaw_bool, bamp, max_step_mpc, cell_len_mpc, delta_einj_list, set_cosmo, max_traj_mpc, use_bp, iter_num=None, tol=None):
    """Write and make output directory which stores relevanat simulation info. The base output direction is simply 'output' in the current working directory.
    Note: the output type (e.g. trajectory) is saved in the filename.

    Parameters
    ----------
    Defined in `simcrpropa_main` docstring.
    """
    # Dictionaries for neat dir names
    pid_dict = {11: "electrons", -11: "positrons", 22: "photons"}
    powlaw_dict = {True: "powlaw", False: "delta"}
    cosmo_dict = {True: "set_cosmo", False: "default_cosmo"}
    bp_dict = {True: "bp", False: "ck"}
    if len(delta_einj_list) == 1:
        delta_einj = delta_einj_list[0]
        odir = os.path.join("output", 
                            pid_dict[pid], 
                            "break_energy",
                            "bx",
                            f"z{redshift}", 
                            f"inj{num_inj}",
                            powlaw_dict[powlaw_bool],
                            f"{delta_einj/1e12}TeV",
                            f"bamp{bamp}G",
                            bp_dict[use_bp],
                            f"tol{tol}",
                            f"iter{iter_num}",
                            cosmo_dict[set_cosmo],
                            f"max_traj_{max_traj_mpc}Mpc",
                            f"max_step_{max_step_mpc}Mpc",
                            f"grid_{(cell_len_mpc):.1e}Mpc",
                            )
    else:       
        string_list = [str(f"{item/1e12}") for item in delta_einj_list]
        odir = os.path.join("output", 
                            pid_dict[pid], 
                            "break_energy",
                            f"z{redshift}", 
                            f"inj{num_inj}",
                            powlaw_dict[powlaw_bool],
                            "elist", # list of energies ... TODO? use `string_list`?
                            f"max_step_{max_step_mpc}Mpc",
                            f"max_traj_{max_traj_mpc}Mpc",
                            f"grid_{(cell_len_mpc):.1e}Mpc",
                            cosmo_dict[set_cosmo]
                            )
    os.makedirs(odir, exist_ok=True)
    print(odir)

    return odir


def setup_sim_modules_and_output(bfield, 
                         max_step_mpc, 
                         max_traj_mpc, 
                         sim_trajectory3d, 
                         sim_event3d, 
                         odir,
                         use_bp,
                         tol=1e-9,
                         sim_detect_all_event3d=False,
                         dsrc=None,
                         origin=ORIGIN,
                         min_step_pc=0.3):
    """Setup entire simulation module, which includes which outputs to save.

    Parameters
    ----------
    bfield : CRPropa magnetic field class
    odir : str
        Base output directory for simulation file(s)
    dsrc : float
        Comoving distance to source in meters. Radius of observer sphere (used if `sim_event3d`=True)
    origin : CRPropa Vector 3d
        Center of observer sphere (used if `sim_event3d`=True)
    min_step_pc : float
        Min allowed propagation step size in *pc*

    Remainder defined in `simcrpropa_main` docstring.
    """
    sim = ModuleList()

    # "The step-size of the particles in our integrator is s_step = min(rg /5, lmax /20) in order to sufficiently resolve the gyrations and the turbulent fluctuations."
    # ref: Anisotropic cosmic ray diffusion in isotropic Kolmogorov turbulence, Reichherzer et al 2022
    # field, tolerance, min step, max step
    if use_bp:
        sim.add(PropagationBP(bfield, tol, min_step_pc*pc, max_step_mpc*Mpc))
    else:
        sim.add(PropagationCK(bfield, tol, min_step_pc*pc, max_step_mpc*Mpc))
    print(f"Max step: {max_step_mpc} Mpc ({max_step_mpc*Mpc} m)")
    sim.add(FutureRedshift())
    sim.add(SynchrotronRadiation(bfield, True, THINNING))
    # Background fields
    cmb = CMB()
    ebl = IRB_Saldana21()
    for bkg in [cmb, ebl]:
        sim.add(EMInverseComptonScattering(bkg, True, THINNING))
        sim.add(EMPairProduction(bkg, True, THINNING))
        sim.add(EMDoublePairProduction(bkg, True, THINNING))
        sim.add(EMTripletPairProduction(bkg, True, THINNING))
    sim.add(MinimumEnergy(100 * MeV))

    if sim_trajectory3d:
        # Bit beyond where observer sphere would be
        sim.add(MaximumTrajectoryLength(max_traj_mpc*Mpc))
        output_traj = TextOutput(os.path.join(odir, 'trajectory3d.txt'), Output.Trajectory3D)
        output_traj.setEnergyScale(eV)
        output_traj.enable(output_traj.CandidateTagColumn)
        output_traj.enable(output_traj.SerialNumberColumn)
        sim.add(output_traj)
        print(output_traj)
    else:
        # Won't be used
        output_traj = None

    if sim_event3d:
        obs = Observer()
        obs.add(ObserverSurface(Sphere(origin, dsrc)))
        output_3d = TextOutput(os.path.join(odir, 'event3d.txt'), Output.Event3D)
        output_3d.setEnergyScale(eV)
        output_3d.enable(output_3d.WeightColumn) # this is required if thinning > 0
        output_3d.enable(output_3d.CandidateTagColumn) # not needed in this analysis
        output_3d.enable(output_3d.SerialNumberColumn)
        obs.onDetection(output_3d)
        obs.setDeactivateOnDetection(True)
        sim.add(obs)
    else:
        output_3d = None

    if sim_detect_all_event3d:
        obs_all = Observer()
        obs_all.add(ObserverDetectAll())
        output_all = TextOutput(os.path.join(odir, 'events_3d_detect_all.txt'), Output.Event3D)
        output_all.enable(output_all.CandidateTagColumn)
        output_all.enable(output_all.SerialNumberColumn)
        obs_all.onDetection(output_all)
        sim.add(obs_all)
    else:
        output_all = None

    sim.setShowProgress(True)
    
    return sim, output_traj, output_3d, output_all


def simcrpropa_main(redshift, 
                    max_traj_mpc,
                    pid, 
                    max_step_mpc, 
                    num_inj, 
                    sim_trajectory3d,  
                    bamp, 
                    delta_einj_list,
                    iter_num,
                    tol,
                    axlims_mpc,
                    use_bp,
                    powlaw_bool=False, 
                    sim_event3d=False, 
                    sim_detect_all_event3d=False,
                    cell_len_mpc=50,
                    set_cosmo=True):
    """Do CRPropa simulation by running all functions in this script.

    Parameters
    ----------
    redshift : float
        Inject particle at this redshift
    max_traj_mpc : float
        Max trajectory distance to track in Mpc
    pid : int
        Particle ID
    max_step : float
        Max propagation step size in Mpc
    num_inj : list[int]
        Number of particles of type `pid` to inject. This can be a one element list.
    sim_trajectory_3d : bool
        If True, save trajectory of particles in simulation using CRPropa's Output.Trajectory3D
    bamp : float
        Fixed amplitude of magnetic field vectors at grid points in Gauss
    delta_einj_list : list[float]
        List of particle energies to inject in eV. Could be a one element list. 
        It *should* match dimensions of `num_inj`? but this is not checked. When setting up the source, `num_inj` is not used, and is not used until `sim.run()`
    iter_num : int
        Iterator number, if running multiple simulations by multiple calls to this function in a loop.
    tol : float
        Tolerance for propagation module. In kpc?
    use_bp : bool
        If True, use BP propagation module.
        If False, use CK propagation module.
    powlaw_bool : bool
        If True, injects a powerlaw for energies of particles of type `pid`.
        If False, injects a delta function for particle energies according to `delta_einj_list`
    sim_event3d : bool
            If True, save events which hit observer sphere, using CRPropa's Output.Event3D
    sim_detect_all_event3d : bool
        If True, detect all particles using CRPropa's ObserverDetectAll with Output.Event3D... did not have much luck with this for some reason.. all empty files if detect all is used ...?
    cell_len : float
        Size of each side of base grid cell in Mpc.
    set_cosmo : bool
        If True, apply cosmology of void catalogs.
        If False, use default CRPropa cosmology (h=0.673, Omega_m=0.315).
    """
    if set_cosmo:
        # https://crpropa.github.io/CRPropa3/api/function_group__PhysicsDefinitions_1ga5c42cc1101552a108ece94569e2c42c4.html
        # Void catalog
        setCosmologyParameters(0.71, 0.27)
        # Default check
        # setCosmologyParameters(0.673, 0.315)

     # Source distance in meters calculated by CRPropa
    dsrc = redshift2ComovingDistance(redshift)
    print(f"Source at a distance: {dsrc/Mpc} Mpc")

    odir = get_outdir(pid=pid, num_inj=num_inj, redshift=z, powlaw_bool=powlaw_bool, bamp=bamp, max_step_mpc=max_step_mpc, cell_len_mpc=cell_len_mpc, delta_einj_list=delta_einj_list, 
                      set_cosmo=set_cosmo, max_traj_mpc=max_traj_mpc, use_bp=use_bp, iter_num=iter_num, tol=tol)

    bfield = create_bfield(grid_len_meters=dsrc, bamp=bamp, cell_len_mpc=cell_len_mpc, outdir=odir, axlims_mpc=axlims_mpc)

    source = create_source(redshift, pid, powlaw_bool, delta_einj_list)

    sim, output_traj, output, output_all = setup_sim_modules_and_output(bfield, 
                                                                        max_step_mpc, 
                                                                        max_traj_mpc, 
                                                                        sim_trajectory3d, 
                                                                        sim_event3d,
                                                                        odir=odir,
                                                                        tol=tol,
                                                                        use_bp=use_bp)

    sim.run(source, num_inj, True)

    # Remember to close files
    if sim_event3d:
        output_3d.close()
    if sim_trajectory3d:
        output_traj.close()
    if sim_detect_all_event3d:
        output_all.close()

    return None


def plotBField(bfield, axlims_mpc, outdir):
    """Plot CRPropa magnetic field `bfield` at fixed points (`xyz_mpc` in the function), and limit axes to `axlims_mpc`.
    Save plot in `outdir`.

    Parameters
    ----------
    bfield : CRPropa magnetic field class
    axlims_mpc : float
        Limit axes to a cube of size `axlims_mpc` * `axlims_mpc` * `axlims_mpc` (in Mpc^3) 
    """
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Copy and pasted numbers from Python np.random.seed because this is not thread safe,
    # so different Cartesian positions were being plotted on different plots. This is not the desired behavior.
    # Cartesian coordinates, each list in the array is [x, y, z] coordinates in Mpc
    xyz_mpc = np.array([[1.95254016e-06, 8.60757465e-06, 4.11053504e-06],
                        [ 1.79532732e-06, -3.05380803e-06,  5.83576452e-06],
                        [-2.49651155e-06,  1.56709200e-05,  1.85465104e-05],
                        [-4.66233925e-06,  1.16690015e-05,  1.15579679e-06],
                        [ 2.72178244e-06,  1.70238655e-05, -1.71585577e-05],
                        [-1.65148280e-05, -1.91912641e-05,  1.33047938e-05],
                        [1.11262700e-05, 1.48004859e-05, 1.91447337e-05],
                        [ 1.19663426e-05, -1.54082551e-06,  1.12211671e-05],
                        [-1.52690230e-05,  5.59684085e-06, -1.42658685e-05],
                        [ 1.77867567e-05,  8.73932870e-07, -3.41352240e-06],
                        [-9.41777552e-06,  1.09693476e-05, -1.75398671e-06],
                        [ 2.73735795e-06, -1.92484080e-05,  4.70541988e-06],
                        [4.48382891e-06, 4.67735987e-06, 1.77499231e-05],
                        [ 7.27281196e-06, -5.61968398e-06, -2.51872185e-06],
                        [ 7.90524784e-06, -1.75909811e-05,  6.67066862e-06],
                        [ 6.82551478e-06, -1.15846976e-05, -1.48429481e-05],
                        [-7.38286596e-06, -5.45156916e-06,  2.80787082e-06],
                        [-2.45593946e-06,  1.95349535e-05, -1.59182076e-05],
                        [-1.16449298e-05, -1.35476193e-05,  6.12433302e-06],
                        [-9.86833590e-06, -1.34756909e-06, -1.02229763e-05]])

    # Point p: (x, y, z)
    for p in xyz_mpc:
        # Convert implied Mpc to meters (CRPropa uses SI units)
        pMeters = Vector3d(p[0], p[1], p[2]) * Mpc 
        field_at_p = bfield.getField(pMeters)
        # Plot in Mpc
        q = ax.quiver(pMeters.x/Mpc, pMeters.y/Mpc, pMeters.z/Mpc, 
                      field_at_p.x, field_at_p.y, field_at_p.z, 
                      normalize=True, length=0.00001)
    ax.set_box_aspect([1, 1, 1])
    lp = 20
    ax.set_xlabel('\nComoving x [Mpc]', labelpad=lp)
    ax.set_ylabel('\nComoving y [Mpc]', labelpad=lp)
    ax.set_zlabel('\nComoving z [Mpc]', labelpad=lp)
    ax.set_xlim([-axlims_mpc, axlims_mpc])
    ax.set_ylim([-axlims_mpc, axlims_mpc])
    ax.set_zlim([-axlims_mpc, axlims_mpc])
    plt.savefig(os.path.join(outdir, 'bfield.png'))
    print(f"--->{os.path.join(outdir, 'bfield.png')}")
    plt.close(fig)

    return None


def initFixedDirectionField(vgrid, Bamplitude):
    """Fill magnetic field grid `vgrid` with fixed direction field (B, B, B) and fixed amplitude `Bamplitude` in Gauss.
    This code based on code from simCRpropa: https://github.com/me-manu/simCRpropa/blob/b3f39b5c77c6b97d19f7db387427d857690444d2/simCRpropa/sim_crpropa.py
    """
    gridArray = vgrid.getGrid()
    nx = vgrid.getNx()
    ny = vgrid.getNy()
    nz = vgrid.getNz()
    print("vgrid: nx = {0:n}, ny = {0:n}, nz = {0:n}".format(
        nx,ny,nz))
    for xi in range(0,nx):
        for yi in range(0,ny):
            for zi in range(0,nz):
                vect3d = vgrid.get(xi,yi,zi)
                # The direction of the source cone is +x, +y
                x = 1
                y = 0
                z = 0
                d = np.sqrt(x*x+y*y+z*z)
                # Unit vectors
                vect3d.x = Bamplitude * x/d
                vect3d.y = Bamplitude * y/d
                vect3d.z = Bamplitude * z/d
    return None


def initRandomField(vgrid, Bamplitude, seed=2308):
    """Fill magnetic field grid `vgrid` with random directions and fixed amplitude `Bamplitude` in Gauss.
    This code is from simCRpropa: https://github.com/me-manu/simCRpropa/blob/b3f39b5c77c6b97d19f7db387427d857690444d2/simCRpropa/sim_crpropa.py
    """
    # Tried two random seeds; both not thread-safe.
    # rng = np.random.default_rng(seed)
    np.random.seed(seed)
    gridArray = vgrid.getGrid()
    nx = vgrid.getNx()
    ny = vgrid.getNy()
    nz = vgrid.getNz()
    print("vgrid: nx = {0:n}, ny = {0:n}, nz = {0:n}".format(
        nx,ny,nz))
    for xi in range(0,nx):
        for yi in range(0,ny):
            for zi in range(0,nz):
                vect3d = vgrid.get(xi,yi,zi)
                x = np.random.uniform(-1,1)
                y = np.random.uniform(-1,1)
                z = np.random.uniform(-1,1)
                # x = rng.uniform(-1, 1)
                # y = rng.uniform(-1, 1)
                # z = rng.uniform(-1, 1)
                d = np.sqrt(x*x+y*y+z*z)
                vect3d.x = Bamplitude * x/d
                vect3d.y = Bamplitude * y/d
                vect3d.z = Bamplitude * z/d
    return None


if __name__ == "__main__":
    # Iterator, to check impact of np.random.seed when using multiple threads
    for itr in [0]:
        # Test difference max step sizes in propagation module
        for i, max_step in enumerate(MAX_STEP_MPC):
            for max_traj in MAX_TRAJ_MPC:
                for z in REDSHIFTS:
                    # Use BP (True), and then CK (False)
                    for bp in [True]:
                        simcrpropa_main(redshift=z, 
                                        pid=PID, 
                                        max_step_mpc=max_step, 
                                        max_traj_mpc=max_traj,
                                        num_inj=NINJ, 
                                        sim_trajectory3d=TRAJECTORY,  
                                        powlaw_bool=False, 
                                        bamp=BAMP, 
                                        delta_einj_list=[E_INJ],
                                        tol=TOL,
                                        axlims_mpc=B_AXLIMS,
                                        iter_num=itr,
                                        use_bp=bp
                                        )

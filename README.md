# `simcrpropa_step_size`

Explore impact of maximum step size on electron trajectory in magnetic field $B = 10^{-13}$ G where $\hat{B} = \dfrac{\hat{x} + \hat{y} + \hat{z}}{\sqrt{3}}$, in order to ensure the maximum step size is set appropriately such that it can resolve the Larmor radius of electrons and positrons.


Tested step sizes: 50, 10, 5, 1, 0.1, $10^{-6}$ Mpc. All resolve the Larmor radius for the lower energy electrons in our simulations (the smallest tracked energy is 100 MeV).

The start of the trajectory of the electron for various max step sizes:
![Beginning part of the trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. The trajectory is confined to a cube centered at (0, 0, 0) with each side of length 2e-5 Mpc. One trajectory is shown per max step size. The maximum tracked trajectory was 0.001 Mpc. The trajectories overlap.](output/trajectory_compare/z0.1/traj_0.001Mpc/stack/bp/start_trajectory3d_prim.png)

The end of the trajectory of the electron for various max step sizes (the Larmor radius is resolved and the trajectories are out of phase):
![Last part of the trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. The last 200 steps of the trajectory are shown. One trajectory is shown per max step size. The maximum tracked trajectory was 0.001 Mpc. The trajectories overlap.](output/trajectory_compare/z0.1/traj_0.001Mpc/stack/bp/end_trajectory3d_prim_n200.png)

The full trajectory (the maximum tracked trajectory is 1 Mpc):
![The full 1 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown per max step size. The trajectories are overlapping straight lines because on these scales the gyromotion cannot be resolved](output/trajectory_compare/z0.1/traj_1Mpc/stack/bp/full_trajectory3d_prim.png)

Near the end of the trajectory, where $r = \sqrt{x^2 + y^2 + z^2}$ the Larmor radius is somewhat resolve, and the trajectories are out of phase:
![Snippet of 1 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown per max step size.](output/trajectory_compare/z0.1/traj_1Mpc/stack/bp/trajectory3d_prim_0.5_to_0.5003.png)

Just one of the lines on the plot above is shown below:
![Snippet of 1 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown, for max step size of 0.1 Mpc.](output/electrons/z0.1/inj1/delta/0.00011TeV/bamp1e-13G/bp/tol1e-09/iter0/set_cosmo/max_traj_1Mpc/max_step_0.1Mpc/grid_5.0e+01Mpc/trajectory3d_prim_0.5_to_0.5003.png)


Tracking a trajectory of 50 Mpc (which is the grid cell size):
Max step size 1 Mpc:
![Last 10000 steps of 50 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown, for max step size of 1 Mpc.](output/electrons/z0.1/inj1/delta/0.00011TeV/bamp1e-13G/bp/tol1e-09/iter0/set_cosmo/max_traj_50Mpc/max_step_1Mpc/grid_5.0e+01Mpc/end_trajectory3d_prim_n10000.png)
![Last 200 steps of 50 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown, for max step size of 1 Mpc.](output/electrons/z0.1/inj1/delta/0.00011TeV/bamp1e-13G/bp/tol1e-09/iter0/set_cosmo/max_traj_50Mpc/max_step_1Mpc/grid_5.0e+01Mpc/end_trajectory3d_prim_n200.png)

Max step size 50 Mpc:
![Last 10000 steps of 50 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown, for max step size of 50 Mpc.](output/electrons/z0.1/inj1/delta/0.00011TeV/bamp1e-13G/bp/tol1e-09/iter0/set_cosmo/max_traj_50Mpc/max_step_50Mpc/grid_5.0e+01Mpc/end_trajectory3d_prim_n10000.png)
![Last 200 steps of 50 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown, for max step size of 50 Mpc.](output/electrons/z0.1/inj1/delta/0.00011TeV/bamp1e-13G/bp/tol1e-09/iter0/set_cosmo/max_traj_50Mpc/max_step_50Mpc/grid_5.0e+01Mpc/end_trajectory3d_prim_n200.png)
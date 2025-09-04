# `simcrpropa_step_size`

## Case: $\hat{B} = B\hat{x}, \vec{v} = c (\hat{x} + \hat{y})$

One electron is in a magnetic field $\hat{B} = B\hat{x}$ with $B = 10^{-13}$ G and the electron has initial velocity $\vec{v} = c (\hat{x} + \hat{y})$.  
The trajectory is tracked up to 1 Mpc.
The Larmor radius is resolved for all tested maximum step sizes. 
The expected Larmor radius of a 110 MeV electron in $B = 10^{-13}$ G is about 1 pc, so a maximum step size of 1 pc is expected to resolve the Larmor radius.

### 110 MeV electron

For 110 MeV electron, the expected Larmor radius is about 1 pc.
It is resolved for all tested maximum step sizes.
The trajectory is tracked up to 1 Mpc.

![Binned trajectory of primary in YX plane for various step sizes for 110 MeV electron](output/trajectory_compare/bx/z0.1/0.00011TeV/traj_1Mpc/trajectory_plots_binned/prim_xy_n100.png)
![Binned trajectory of primary in ZX plane for various step sizes for 110 MeV electron](output/trajectory_compare/bx/z0.1/0.00011TeV/traj_1Mpc/trajectory_plots_binned/prim_zx_n100.png)
![Binned trajectory of primary in ZY plane for various step sizes for 110 MeV electron](output/trajectory_compare/bx/z0.1/0.00011TeV/traj_1Mpc/trajectory_plots_binned/prim_zy_n100.png)

![Energy of primary 110 MeV electron versus comoving X various step sizes.](output/trajectory_compare/bx/z0.1/0.00011TeV/traj_1Mpc/energy_plots/energy_vs_x_prim.png)

Note: 
* The electron at 110 MeV is just above the simulation break energy (100 MeV) so no noticablle evolution of the radius is expected to be observed with this particular simulation, because as the radius decreases so does the energy $E \propto r_L$ at which point the particle may be untracked.
* Many steps are needed to traverse the full perimeter.

### 1 TeV electron 

For 1 TeV electron, the expected Larmor radius is about 0.01 Mpc.
It is resolved for all tested maximum step sizes.
The trajectory is tracked up to 1 Mpc.

![Binned trajectory of primary in YX plane for various step sizes for 1 TeV electron](output/trajectory_compare/bx/z0.1/1.0TeV/traj_1Mpc/trajectory_plots_binned/prim_xy_n200.png)
![Binned trajectory of primary in ZX plane for various step sizes for 1 TeV electron](output/trajectory_compare/bx/z0.1/1.0TeV/traj_1Mpc/trajectory_plots_binned/prim_zx_n200.png)
![Binned trajectory of primary in ZY plane for various step sizes for 1 TeV electron](output/trajectory_compare/bx/z0.1/1.0TeV/traj_1Mpc/trajectory_plots_binned/prim_zy_n200.png)

![Energy of primary 1 TeV electron versus comoving X various step sizes.](output/trajectory_compare/bx/z0.1/1.0TeV/traj_1Mpc/energy_plots/energy_vs_x_prim.png)

Note:
* Secondaries are expected to be tracked in this case because they do not quickly fall below the break energy of 100 MeV. 
* The electron loses energy via interactions so its energy and therefore Larmor radius are not fixed.
* The creation of cascades is via Monte Carlo methods so results are not expected to be identical.

---

### Cascades

### 110 MeV electron

![Binned trajectory of cascade in YX plane for various step sizes for 110 MeV electron](output/trajectory_compare/bx/z0.1/0.00011TeV/traj_1Mpc/trajectory_plots_binned/casc_xy_n100.png)
![Binned trajectory of cascade in ZX plane for various step sizes for 110 MeV electron](output/trajectory_compare/bx/z0.1/0.00011TeV/traj_1Mpc/trajectory_plots_binned/casc_zx_n100.png)
![Binned trajectory of cascade in ZY plane for various step sizes for 110 MeV electron](output/trajectory_compare/bx/z0.1/0.00011TeV/traj_1Mpc/trajectory_plots_binned/casc_zy_n100.png)

![Energy of secondary particles versus comoving X various step sizes.](output/trajectory_compare/bx/z0.1/0.00011TeV/traj_1Mpc/energy_plots/energy_vs_x_casc.png)

### 1 TeV electron 

![Binned trajectory of cascade in YX plane for various step sizes for 1 TeV electron](output/trajectory_compare/bx/z0.1/1.0TeV/traj_1Mpc/trajectory_plots_binned/casc_xy_n200.png)
![Binned trajectory of cascade in ZX plane for various step sizes for 1 TeV electron](output/trajectory_compare/bx/z0.1/1.0TeV/traj_1Mpc/trajectory_plots_binned/casc_zx_n200.png)
![Binned trajectory of cascade in ZY plane for various step sizes for 1 TeV electron](output/trajectory_compare/bx/z0.1/1.0TeV/traj_1Mpc/trajectory_plots_binned/casc_zy_n200.png)

![Energy of secondary particles versus comoving X various step sizes.](output/trajectory_compare/bx/z0.1/1.0TeV/traj_1Mpc/energy_plots/energy_vs_x_casc.png)

---

### Case: $\hat{B} = \dfrac{\hat{x} + \hat{y} + \hat{z}}{\sqrt{3}}, \vec{v} = c \hat{x}$

Explore impact of maximum step size on electron trajectory in magnetic field $B = 10^{-13}$ G where $\hat{B} = \dfrac{\hat{x} + \hat{y} + \hat{z}}{\sqrt{3}}$, in order to ensure the maximum step size is set appropriately such that it can resolve the Larmor radius of electrons and positrons.


Tested step sizes: 50, 10, 5, 1, 0.1, $10^{-6}$ Mpc. All resolve the Larmor radius for the lower energy electrons in our simulations (the smallest tracked energy is 100 MeV).

The start of the trajectory of the electron for various max step sizes:
![Beginning part of the trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. The trajectory is confined to a cube centered at (0, 0, 0) with each side of length 2e-5 Mpc. One trajectory is shown per max step size. The maximum tracked trajectory was 0.001 Mpc. The trajectories overlap.](output/trajectory_compare/bx_by_bz/z0.1/traj_0.001Mpc/stack/bp/start_trajectory3d_prim.png)

The end of the trajectory of the electron for various max step sizes (the Larmor radius is resolved and the trajectories are out of phase):
![Last part of the trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. The last 200 steps of the trajectory are shown. One trajectory is shown per max step size. The maximum tracked trajectory was 0.001 Mpc. The trajectories overlap.](output/trajectory_compare/bx_by_bz/z0.1/traj_0.001Mpc/stack/bp/end_trajectory3d_prim_n200.png)

The full trajectory (the maximum tracked trajectory is 1 Mpc):
![The full 1 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown per max step size. The trajectories are overlapping straight lines because on these scales the gyromotion cannot be resolved](output/trajectory_compare/bx_by_bz/z0.1/traj_1Mpc/stack/bp/full_trajectory3d_prim.png)

Near the end of the trajectory, where $r = \sqrt{x^2 + y^2 + z^2}$ the Larmor radius is somewhat resolve, and the trajectories are out of phase:
![Snippet of 1 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown per max step size.](output/trajectory_compare/bx_by_bz/z0.1/traj_1Mpc/stack/bp/trajectory3d_prim_0.5_to_0.5003.png)

Just one of the lines on the plot above is shown below:
![Snippet of 1 Mpc trajectory of a 110 MeV electron at redshift 0.1 in the magnetic field. One trajectory is shown, for max step size of 0.1 Mpc.](output/electrons/break_energy/bx_by_bz/z0.1/inj1/delta/0.00011TeV/bamp1e-13G/bp/tol1e-09/iter0/set_cosmo/max_traj_1Mpc/max_step_0.1Mpc/grid_5.0e+01Mpc/trajectory3d_prim_0.5_to_0.5003.png)

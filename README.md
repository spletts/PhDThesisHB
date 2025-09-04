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
"""Utilities for plotting (primary) particle trajectories (ignoring cascade particles)"""


import matplotlib.pyplot as plt
import os
from numpy import genfromtxt, unique
import numpy as np
from crpropa import *
from collections import Counter
from astropy import units as u
from astropy.table import Table
import sim


# Particle IDs
PID_DICT = {11: "electron", -11: "positron", 22: "photon"}
# Create list (rather than iterator) of markers and colors so each file in list of filenames is indexed consistently across multiple plots
MARKERS = ['o', 's', 'd', 'h', '>', '^', 'X']
TAB10_COLORS = [
    "#1f77b4",  # Blue
    "#ff7f0e",  # Orange
    "#2ca02c",  # Green
    "#d62728",  # Red
    "#9467bd",  # Purple
    "#8c564b",  # Brown
    "#e377c2",  # Pink
    "#7f7f7f",  # Gray
    "#bcbd22",  # Olive
    "#17becf"   # Cyan
    ]
# For getting output filename. Must match `bp_dict` in the function `sim.get_outdir`
BP_DICT = {True: "bp", False: "ck"}
# For plot display
BP_PLOT_DICT = {True: "Borish Push", False: "Cash-Karp"}


def xyz_last_slice_of_large_file(fn, num=200):
    """Efficiently read last `num` lines of large file `fn` *which correspond to primary particles*, and return the x, y, z comoving coordinates in Mpc.
    Efficiently read file by reading it in reverse order with `tac`, rather than reading with Python and storing the entire file in memory.

    Parameters
    ----------
    fn : str 
        CRPropa Output.Trajectory3D
        Header of `fn`: #	D	SN	ID	E	X	Y	Z	Px	Py	Pz	SN0	SN1	tag
    num : int
        Read the last `num` lines of the file `fn` *which contain primary particles*

    Returns
    -------
    x, y, z: array_like[float]
        Comoving coordinates in Mpc
    """
    import subprocess

    # Indices/columns for X, Y, and Z in header 
    ix, iy, iz = 4, 5, 6
    _x, _y, _z = [], [], []
    # Primary particles tagged as PRIM. 
    # `tac` to reverse search. Add `| tac`` at end to put back in original order
    cmd = f'tac "{fn}" | grep -m {num} "PRIM" | tac'
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True).stdout
    # Get each line as element in list
    result_lines = result.splitlines()
    for rl in result_lines:
        # Split each line into list, with each element a cell in that row
        rl_list = rl.split("\t")       
        _x.append(float(rl_list[ix]))
        _y.append(float(rl_list[iy]))
        _z.append(float(rl_list[iz]))
    x, y, z = np.array(_x), np.array(_y), np.array(_z)

    return x, y, z


def xyz_by_distance_of_large_file(fn, r1, r2, max_num=200, comment_char="#"):
    """Efficiently read large file `fn` line by line so it is not stored in memory, and return the relevant x, y, z comoving coordinates in Mpc which satisfy a distance cut. 
    Specifically, at *most*, return the first `max_num` *primary* particle coordinates (comoving, in Mpc) that satisfy:
    r1 <= (x^2 + y^2 + z^2)^(1/2) <= r2 (Mpc), so the stacked plot is over the same axes range. Cascade particles are ignored.

    Efficiently read file by reading it line, rather than reading with Python and storing the entire file in memory. 
    For very large files, the job is OOM killed when requesting >60GB of memory. It is also slow.
    This function is similar to `last_xyz_of_large_file`, but I found that when slicing the last `num` lines of different files, they were at different 
    enough x, y, z coordinates (why? the same maximum tracked trajectory was set) that the spiral/gyroradius could not be resolved due to the x, y, z extent of the plot.
    So this function looks for all points within a specific distance range, where the distance is r^2 = x^2 + y^2 + z^2. For the assumed magnetic field (B, B, B), the particle spirals 
    around this line for all steps, so distance is calculated as x^2 + y^2 + z^2 because the particle is in a magnetic field with fixed direction.
    
    Parameters
    ----------
    fn : str 
        CRPropa Output.Trajectory3D
        Header of `fn`: #	D	SN	ID	E	X	Y	Z	Px	Py	Pz	SN0	SN1	tag
    r1 : float
        Lower bound on distance sqrt(x^2 + y^2 + z^2) in Mpc
    r2 : float
        Upper bound on distance sqrt(x^2 + y^2 + z^2) in Mpc.
        Note that a good value for r2 - r1 is ~10^-5 Mpc
    max_num : int
        At *most* return the first `max_num` *primary* coordinates between `r1` and `r2`. 
        If there are less than `max_num` points between `r1` and `r2`, then they will all be returned.
    comment_char : str
        Lines that start with `comment_char` in `fn` are comments and will be ignored when reading the file.

    Returns
    -------
    x, y, z: array_like[float]
        Comoving coordinates in Mpc
    """
    # Indices/columns for X, Y, and Z in header 
    ix, iy, iz = 4, 5, 6
    _x, _y, _z = [], [], []
    # Number of coordinates in the lists `_x`, `_y` and `_z`
    coord_counter = 0
    with open(fn, 'r') as f:
        for line in f:
            # Ensure line corresponds to a primary particle, and is not a comment line
            if 'PRIM' in line and line.startswith(comment_char) is False and coord_counter <= max_num:
                l = line.split("\t")
                x, y, z = float(l[ix]), float(l[iy]), float(l[iz])
                r = np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))
                if r <= r2 and r >= r1 and coord_counter <= max_num:
                    coord_counter += 1
                    _x.append(x)
                    _y.append(y)
                    _z.append(z)
    x, y, z = np.array(_x), np.array(_y), np.array(_z)

    return x, y, z


def stack_full_3d_trajectories(fn_data_list, max_step_list, e_inj, redshift, traj, use_bp, ii=None, n_inj="one"):
    """Plot all 3d trajectories in `fn_data_list` on one plot. 
    Plot the full span of trajectory, or if `ii` is not None, plot ever `ii` step in the trajectory.
    
    Parameters
    ----------
    fn_data_list : list[str]
        This can be a one element list to plot one trajectory.
        Each element is a CRPropa Output.Trajectory3D
        with header: #	D	SN	ID	E	X	Y	Z	Px	Py	Pz	SN0	SN1	tag
    max_step_list : str[float]
        This can be a one element list to plot one trajectory.
        Max step allowed in Mpc for each element in `fn_data_list`. 
        The first element in `fn_data_list` corresponds to the first element in `max_step_list`.
    redshift : float
    e_inj : float
        Injected energy in eV
    traj : float or str
        Max tracked trajectory. Used for output directory name (so can be float or str).
    use_bp : bool
        If True, BP propagation was used in all simulations in `fn_data_list`.
        IF False, CK method was used.
    ii : int or None
        If not None, plot every `ii` point.
    n_inj : str of float
        Number of injected particles. This is used for display so can be a str or float.
    """

    # Setup axes: one for primary
    # The cascade particle trajectories are random in the sense that whether an interaction occurs (and what interaction occurs) 
    # at a specific step is determined using 2 random numbers
    fig_prim = plt.figure(figsize=(14, 11))
    ax_prim = fig_prim.add_axes(111, projection="3d")
    # Still keeping infrastructure for cascade
    axes = [ax_prim]
    figs = [fig_prim]

    if len(fn_data_list) != 1:
        odir = os.path.join("output", "trajectory_compare", f"z{redshift}", f"traj_{traj}Mpc", "stack", BP_DICT[use_bp])
    else:
        # If plotting one trajectory, save plot in the same directory as the trajectory file `fn_data_list[0]`
        odir = os.path.dirname(fn_data_list[0])

    os.makedirs(odir, exist_ok=True)
    fn_plot_list = [os.path.join(odir, "full_trajectory3d_prim.png"),]   

    for k, fn_data in enumerate(fn_data_list):
        data = genfromtxt(fn_data, names=True, dtype=None)
        tags = data["tag"][:]
        pids = data["ID"][:]
        idx_prim = np.where(tags == b'PRIM')[0]
        # Only interested in primary 
        for j, idx in enumerate([idx_prim]):
            # Primary
            ax = axes[j]
            fig = figs[j]
            data_mask = data[idx]
            pid_mask = pids[idx]
            x_arr, y_arr, z_arr = data_mask['X'], data_mask['Y'], data_mask['Z']
            if ii is not None:
                x_arr, y_arr, z_arr  = x_arr[::ii], y_arr[::ii], z_arr[::ii]
            ax.scatter(x_arr, y_arr, z_arr, 
                   label=f"Max step: {max_step_list[k]} Mpc", 
                   facecolor="none", marker=MARKERS[k], color=TAB10_COLORS[k])

    ax.scatter(0, 0, 0, marker='*', color="k", label=f"Injection location", zorder=100)
    # Make room for legend
    fig.subplots_adjust(right=0.8)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    lp = 25
    ax.set_xlabel('Comoving x [Mpc]', labelpad=lp)
    ax.set_ylabel('Comoving y [Mpc]', labelpad=lp)
    ax.set_zlabel('Comoving z [Mpc]', labelpad=lp)
    ax.set_aspect('equal')
    if ii is not None:
        ax.set_title(rf"Trajectory (every {ii} steps) of {n_inj} {e_inj/1e6} MeV {PID_DICT[pid_mask[0]]} at $z=${redshift} in $B=10^{{-13}}$ G grid" + f"\n(propagation: {BP_PLOT_DICT[use_bp]})")
    else:
        ax.set_title(rf"Trajectory of {n_inj} {e_inj/1e6} MeV {PID_DICT[pid_mask[0]]} at $z=${redshift} in $B=10^{{-13}}$ G grid" + f"\n(propagation: {BP_PLOT_DICT[use_bp]})")
    for j, fig in enumerate(figs):
        fig.savefig(fn_plot_list[j])
        print(f"Plotted: {fn_plot_list[j]}")
        plt.close(fig)

    return None


def stack_last_slice_3d_trajectories(fn_data_list, max_step_list, e_inj, redshift, traj, num, use_bp, n_inj="one"):
    """Plot the last `num` steps of all 3d trajectories in `fn_data_list` on one plot. 
    
    Parameters
    ----------
    See `stack_full_3d_trajectories` docstring.
    """

    fig_prim = plt.figure(figsize=(14, 11))
    ax_prim = fig_prim.add_axes(111, projection="3d")
    axes = [ax_prim]
    figs = [fig_prim]

    if len(fn_data_list) != 1:
        odir = os.path.join("output", "trajectory_compare", f"z{redshift}", f"traj_{traj}Mpc", "stack", BP_DICT[use_bp])
    else:
        # If plotting one trajectory, save plot in the same directory as the trajectory file `fn_data_list[0]`
        odir = os.path.dirname(fn_data_list[0])
    os.makedirs(odir, exist_ok=True)
    fn_plot_list = [os.path.join(odir, f"end_trajectory3d_prim_n{num}.png"),]   

    # Primary only
    ax = axes[0]
    fig = figs[0]
    for k, fn_data in enumerate(fn_data_list):
        # This function looks at primary particles only
        x_arr, y_arr, z_arr = xyz_last_slice_of_large_file(fn_data, num)
        ax.scatter(x_arr, y_arr, z_arr, 
                   label=f"Max step: {max_step_list[k]} Mpc", 
                   facecolor="none", marker=MARKERS[k], color=TAB10_COLORS[k])
    # Make room for legend
    fig.subplots_adjust(right=0.8)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    lp = 25
    ax.set_xlabel('Comoving x [Mpc]', labelpad=lp)
    ax.set_ylabel('Comoving y [Mpc]', labelpad=lp)
    ax.set_zlabel('Comoving z [Mpc]', labelpad=lp)
    ax.set_aspect('equal')
    ax.set_title(rf"Trajectory (last {num} steps of {traj} Mpc) of {n_inj} {e_inj/1e6} MeV electron at $z=${redshift} in $B=10^{{-13}}$ G grid" +  f"\n(propagation: {BP_PLOT_DICT[use_bp]})")
    for j, fig in enumerate(figs):
        fig.savefig(fn_plot_list[j])
        print(f"Plotted: {fn_plot_list[j]}")
        plt.close(fig)

    return None


def stack_slice_by_distance_3d_trajectories(fn_data_list, max_step_list, e_inj, redshift, traj, max_num, use_bp, r1=0.2, r2=0.2002, n_inj="one"):
    """Plot at most `max_num` points between the distance r1 <= r <= r2 in the 3d trajectories of all files in `fn_data_list` on one plot. 
    The distance r^2 = x^2 + y^2 + z^2. See `xyz_by_distance_of_large_file`.
    
    Parameters
    ----------
    See `stack_full_3d_trajectories` and `xyz_by_distance_of_large_file` docstrings.
    """

    fig_prim = plt.figure(figsize=(14, 11))
    ax_prim = fig_prim.add_axes(111, projection="3d")
    axes = [ax_prim]
    figs = [fig_prim]

    if len(fn_data_list) != 1:
        odir = os.path.join("output", "trajectory_compare", f"z{redshift}", f"traj_{traj}Mpc", "stack", BP_DICT[use_bp])
    else:
        odir = os.path.dirname(fn_data_list[0])
    os.makedirs(odir, exist_ok=True)
    fn_plot_list = [os.path.join(odir, f"trajectory3d_prim_{r1}_to_{r2}.png"),]   

    # Primary only
    ax = axes[0]
    fig = figs[0]
    for k, fn_data in enumerate(fn_data_list):
        # This function looks at primary particles only
        x_arr, y_arr, z_arr = xyz_by_distance_of_large_file(fn_data, r1, r2, max_num)
        ax.scatter(x_arr, y_arr, z_arr, 
                   label=f"Max step: {max_step_list[k]} Mpc", 
                   facecolor="none", marker=MARKERS[k], color=TAB10_COLORS[k])
    # Make room for legend
    fig.subplots_adjust(right=0.8)
    ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
    lp = 25
    ax.set_xlabel('Comoving x [Mpc]', labelpad=lp)
    ax.set_ylabel('Comoving y [Mpc]', labelpad=lp)
    ax.set_zlabel('Comoving z [Mpc]', labelpad=lp)
    ax.set_aspect('equal')
    ax.set_title(rf"Trajectory (first {max_num} steps between {r1} < $r$ < {r2} Mpc) of {n_inj} {e_inj/1e6} MeV electron in $B=10^{{-13}}$ G grid")
    ax.xaxis.get_major_formatter().set_useOffset(False)
    ax.yaxis.get_major_formatter().set_useOffset(False)
    ax.zaxis.get_major_formatter().set_useOffset(False)
    for j, fig in enumerate(figs):
        fig.savefig(fn_plot_list[j])
        print(fn_plot_list[j])
        plt.close(fig)

    return None


def stack_confine_3d_trajectory(fn_data_list, e_inj, redshift, max_step_list, traj, use_bp, n_inj="one", axlim=1e-5):
    """Plot beginning of 3d trajectories of all files in `fn_data_list` on one plot, and limit/mask trajectories to a cube within:
    -`axlim` <= x <= `axlim`, -`axlim` <= y <= `axlim`, -`axlim` <= y <= `axlim`, where `axlim` is in Mpc.

    Parameters
    ----------
    axlim : float
        Confine 3d axis to a cube, each of side length 2*`axlim` Mpc, spanning -`axlim` to +`axlim`.
    See `stack_full_3d_trajectories` docstring.
    """

    fig_prim = plt.figure(figsize=(14, 11))
    ax_prim = fig_prim.add_axes(111, projection="3d")
    axes = [ax_prim]
    figs = [fig_prim]
    fn_plot_list = [os.path.join("output", "trajectory_compare", f"z{redshift}", f"traj_{traj}Mpc", "stack", BP_DICT[use_bp], "start_trajectory3d_prim.png"),]  
    os.makedirs(os.path.dirname(fn_plot_list[0]), exist_ok=True)

    for i, fn_data in enumerate(fn_data_list):
        data = genfromtxt(fn_data, names=True, dtype=None)
        pids = data["ID"][:]
        tags = data["tag"][:]
        idx_prim = np.where(tags == b'PRIM')[0]
        for j, idx in enumerate([idx_prim]):
            data_mask = data[idx]
            pid_mask = pids[idx]
            uniq_pid_mask = np.unique(pid_mask)
            ax = axes[j]
            fig = figs[j]
            # Comoving coords in Mpc
            _x_arr, _y_arr, _z_arr = data_mask['X'], data_mask['Y'], data_mask['Z']
            mask_idx = np.where((abs(_x_arr) <= axlim) & (abs(_y_arr) <= axlim) & (abs(_z_arr) <= axlim))[0]
            x_arr, y_arr, z_arr = _x_arr[mask_idx], _y_arr[mask_idx], _z_arr[mask_idx]
            # Handle primarys (one particle type)
            if j == 0:
                # For primary, all `pid_mask` are the same
                ax.scatter(x_arr, y_arr, z_arr, 
                        alpha=0.8, color=TAB10_COLORS[i], marker=MARKERS[i], facecolor="none",
                        label=f"Max step: {max_step_list[i]} Mpc")
            else:
                # Plot full trajectory
                print(f"--->Number of cascade particles: {len(idx)}")
                ax.scatter(_x_arr, _y_arr, _z_arr, label="Cascade particle(s)", alpha=0.8)             

    for k, ax in enumerate(axes):
        ax.scatter(0, 0, 0, marker='*', color="k", label=f"Injection location")
        # Make room for legend
        fig.subplots_adjust(right=0.8)
        ax.legend(loc='center left', bbox_to_anchor=(1.05, 0.5))
        lp = 25
        ax.set_xlabel('\nComoving x [Mpc]', labelpad=lp)
        ax.set_ylabel('\nComoving y [Mpc]', labelpad=lp)
        ax.set_zlabel('\nComoving z [Mpc]', labelpad=lp)
        # Limit primary axis span, not cascade (if plotting cascade at all)
        if k == 0:
            ax.set_xlim([-axlim, axlim])
            ax.set_ylim([-axlim, axlim])
            ax.set_zlim([-axlim, axlim])
        ax.set_aspect('equal')
        ax.set_title(rf"Beginning trajectory of {n_inj} {e_inj/1e6} MeV electron at $z=${z} in $B=10^{{-13}}$ G grid" + f"\n(propagation: {BP_PLOT_DICT[use_bp]})")  # + "\n" + rf"{COSMO_DICT[default_cosmo]}
    for j, fig in enumerate(figs):
        fig.savefig(fn_plot_list[j])
        print(f"Plotted: {fn_plot_list[j]}")
        plt.close(fig)

    return None


if __name__ == "__main__":

    plt.rcParams.update({'font.size': 14, 'axes.labelsize': 14, 'legend.fontsize': 12})

    # Dictionaries for neater plot labels
    COSMO_DICT = {True: "Default cosmology ($h=0.673, \Omega_m=0.315$)", False: "Nondefault cosmology ($h=0.71, \Omega_m=0.27$)"}

    MAX_STEP_MPC = sim.MAX_STEP_MPC
    E_INJ = sim.E_INJ
    MAX_TRAJ_MPC = sim.MAX_TRAJ_MPC

    # Which plots to make
    # Entire trajectory
    FULL_TRAJ = False
    # Stack multiple simulation outputs on one plot
    STACK_TRAJ = True
    # Plot on simulation output per plot
    INDIVIDUAL_TRAJ = False


    for bp in [True]:
        for z in sim.REDSHIFTS:
            ii_list = [10, 1000]
            for i, max_traj in enumerate(MAX_TRAJ_MPC):
                # Want trajectories with the source at the same redshift and same max trajectory traced on the same plot
                fn_list_per_traj_and_z = []
                for max_step in MAX_STEP_MPC:
                    fn =  f"output/electrons/z{z}/inj1/delta/0.00011TeV/bamp1e-13G/{BP_DICT[bp]}/tol1e-09/iter0/set_cosmo/max_traj_{max_traj}Mpc/max_step_{max_step}Mpc/grid_5.0e+01Mpc/trajectory3d.txt"
                    fn_list_per_traj_and_z.append(fn)
                    print(f"Analysing {fn}")

                    if INDIVIDUAL_TRAJ:
                        if max_traj == 1:
                            r1 = 0.5
                            r2 =  r1 + 0.0003
                            stack_slice_by_distance_3d_trajectories(fn_data_list=[fn], e_inj=E_INJ, redshift=z, max_step_list=[max_step], use_bp=bp,
                                                                    traj=max_traj, r1=r1, r2=r2, max_num=300)
                if STACK_TRAJ:
                ##############################
                    # Out of max step loop; stack list of files in `fn_list_per_traj_and_z`
                    stack_confine_3d_trajectory(fn_data_list=fn_list_per_traj_and_z, e_inj=E_INJ, redshift=z, max_step_list=MAX_STEP_MPC, traj=max_traj, axlim=1e-5, use_bp=bp)
                    if max_traj == 0.001:
                        stack_last_slice_3d_trajectories(fn_data_list=fn_list_per_traj_and_z, max_step_list=MAX_STEP_MPC, e_inj=E_INJ, redshift=z, traj=max_traj, num=200, use_bp=bp)
                    if max_traj == 1:
                        r1 = 0.5
                        r2 =  r1 + 0.0003
                        stack_slice_by_distance_3d_trajectories(fn_data_list=fn_list_per_traj_and_z, e_inj=E_INJ, redshift=z, max_step_list=MAX_STEP_MPC, use_bp=bp,
                                                                traj=max_traj, r1=r1, r2=r2, max_num=300)

                    if FULL_TRAJ:
                        stack_full_3d_trajectories(fn_data_list=fn_list_per_traj_and_z, max_step_list=MAX_STEP_MPC, e_inj=E_INJ, redshift=z, traj=max_traj, use_bp=bp, ii=ii_list[i])

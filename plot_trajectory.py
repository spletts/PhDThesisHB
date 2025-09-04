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
import subprocess


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
CMAPS = ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
# For getting output filename. Must match `bp_dict` in the function `sim.get_outdir`
BP_DICT = {True: "bp", False: "ck"}
# For plot display
BP_PLOT_DICT = {True: "Boris Push", False: "Cash-Karp"}


def plot_energy_versus_x(fn_data_list, max_step_list, redshift, e_inj, traj, use_bp):
    """Plot the particle energy at each step's X coordinate (comoving, in Mpc). Separate the primary and cascade particles (one plot is made for each and saved separately, not in subplots)
    The choice of orientation of the axes and the choice of plotting X is assuming the magnetic field is along +x, and the source velocity is +x, +y.
    """
    fig_prim = plt.figure(figsize=(14, 10))
    ax_prim = fig_prim.add_axes(111)
    fig_casc = plt.figure(figsize=(14, 10))
    ax_casc = fig_casc.add_axes(111)
    axes = [ax_prim, ax_casc]
    figs = [fig_prim, fig_casc]
    titles = ["Primary", "Cascade"]

    if len(fn_data_list) != 1:
        odir = os.path.join("output", "trajectory_compare", "bx", f"z{redshift}",  f"{e_inj/1e12}TeV", f"traj_{traj}Mpc", "energy_plots")
    else:
        # If plotting one trajectory, save plot in the same directory as the trajectory file `fn_data_list[0]`
        odir = os.path.dirname(fn_data_list[0])
    os.makedirs(odir, exist_ok=True)
    fn_plot_list = [os.path.join(odir, "energy_vs_x_prim.png"),
                    os.path.join(odir, "energy_vs_x_casc.png")]   

    for k, fn_data in enumerate(fn_data_list):
        data = genfromtxt(fn_data, names=True, dtype=None)
        tags = data["tag"][:]
        pids = data["ID"][:]
        idx_prim = np.where(tags == b'PRIM')[0]
        idx_casc = np.where(tags != b'PRIM')[0]
        for j, idx in enumerate([idx_prim, idx_casc]):
            ax = axes[j]
            fig = figs[j]
            data_mask = data[idx]
            x_arr, energy = data_mask['X'], data_mask['E']
            ax.scatter(x_arr, energy, 
                        label=f"Max step: {max_step_list[k]} Mpc", 
                        facecolor="none", marker=MARKERS[k], color=TAB10_COLORS[k])
  
    for j, ax in enumerate(axes):
        if j == 0:
            ax.axhline(e_inj, label="Injected energy", color="k", ls="--")
            print("Plotted injected energy")
        ax.set_yscale("log")
        ax.legend()
        ax.set_xlabel('Comoving x [Mpc]')
        ax.set_ylabel('Energy [eV]')
        ax.set_title(titles[j])
    for j, fig in enumerate(figs):
        fig.savefig(fn_plot_list[j])
        print(f"Plotted: {fn_plot_list[j]}")
        plt.close(fig)

    return None


def subplot_binned_2d_trajectories(fn_data_list, max_step_list, redshift, traj, e_inj, plot_prim=True, n=100):
    """Bin full trajectories of the primary particle (`plot_prim=True`) or cascade particles (`plot_prim=False`) in 3 planes: 
    XY, ZY, ZX; each 2d histogram has `n` bins per dimension.
    The entire trajectory distance is plotted.
    There is one figure per plane: XY, ZY, ZX. For each figure, there is one subplots per element in `fn_data_list` (no overplotting).
    If there are two elements in `fn_data_list`, indicating max steps of 1 Mpc and 10 Mpc, there is one figure of the YX planes
    with two subplots, where subplot has the results of max step 1 Mpc and the other subplot has the results of max step 10 Mpc.
    And so on for the ZY, ZX planes.
    """
    # Handle output
    odir = os.path.join("output", "trajectory_compare", "bx", f"z{redshift}", f"{e_inj/1e12}TeV", f"traj_{traj}Mpc", "trajectory_plots_binned")
    os.makedirs(odir, exist_ok=True)
    plot_prim_dict = {True: "prim", False: "casc"}

    # Setup figures and their styles
    if len(fn_data_list)/2 % 2 == 0:
        fig_xy, (axes_xy) = plt.subplots(2, int(len(fn_data_list)/2), figsize=(15, 10), sharex=True, sharey=True)
        fig_zy, (axes_zy) = plt.subplots(2, int(len(fn_data_list)/2), figsize=(15, 10), sharex=True, sharey=True)
        fig_zx, (axes_zx) = plt.subplots(2, int(len(fn_data_list)/2), figsize=(15, 10), sharex=True, sharey=True)
    else:
        fig_xy, (axes_xy) = plt.subplots(1, len(fn_data_list), figsize=(5*len(fn_data_list), 5), sharex=True, sharey=True)
        fig_zy, (axes_zy) = plt.subplots(1, len(fn_data_list), figsize=(5*len(fn_data_list), 5), sharex=True, sharey=True)
        fig_zx, (axes_zx) = plt.subplots(1, len(fn_data_list), figsize=(5*len(fn_data_list), 5), sharex=True, sharey=True)
    figs = [fig_xy, fig_zy, fig_zx]
    axes = [axes_xy.flatten(), axes_zy.flatten(), axes_zx.flatten()]
    units = "Mpc"
    # Format: [x label, y label]
    coord_labels = [ [f"Comoving X [{units}]", f"Comoving Y [{units}]"],
                     [f"Comoving Y [{units}]", f"Comoving Z [{units}]"],
                     [f"Comoving X [{units}]", f"Comoving Z [{units}]"] ]
    plot_titles = ["XY", "ZY", "ZX"]

    for i, fn_data in enumerate(fn_data_list):  
        print(f"{len(fn_data_list)} files to be plotted on {len(axes_xy)} axes")
        data = genfromtxt(fn_data, names=True, dtype=None)
        tags = data["tag"][:]
        # Mask to keep only primaries/secondaries as indicated by `plot_prim`
        if plot_prim:
            idx_prim = np.where(tags == b'PRIM')[0]
            data_mask = data[idx_prim]
        else:
            idx_casc = np.where(tags != b'PRIM')[0]
            data_mask = data[idx_casc]
        # Comoving coords in Mpc
        x_arr, y_arr, z_arr = data_mask['X'], data_mask['Y'], data_mask['Z']
        print(f"Number of coordinates {len(x_arr)}")
        # Format: [x, y]. Get this once per time `fn_data` is read
        coord_pairs = [[x_arr, y_arr],
                       [y_arr, z_arr],
                       [x_arr, z_arr] 
                       ]
        for j, _ in enumerate(['xy', 'zy', 'zx']):
            fig = figs[j]
            ax = axes[j][i]
            h, xedges, yedges = np.histogram2d(coord_pairs[j][0], coord_pairs[j][1], bins=[n, n])
            # Transpose `h` due to this error: TypeError: Dimensions of C (9, 1) should be one smaller than X(10) and Y(2) while using shading='flat' see help(pcolormesh)
            im = ax.pcolormesh(xedges, yedges, np.log10(h.T), label=f"{max_step_list[i]} Mpc")
            fig.colorbar(im, ax=ax, pad=0.01, label="$\log_{10}$(counts)")
            fig.supxlabel(coord_labels[j][0])
            fig.supylabel(coord_labels[j][1])
            ax.set_title(f"Max step: {max_step_list[i]} Mpc")
            # Circle expected in ZY plane
            if j == 1:
                ax.set_aspect('equal')
            fig.suptitle(f"Binned trajectory in {plot_titles[j]} plane")
    planes = ['xy', 'zy', 'zx']
    for k, f in enumerate(figs):
        fn_plot = os.path.join(odir, f"{plot_prim_dict[plot_prim]}_{planes[k]}_n{n}.png")
        f.savefig(fn_plot)
        print(fn_plot)
        plt.close(f)

    return None


if __name__ == "__main__":

    plt.rcParams.update({'font.size': 14, 'axes.labelsize': 14, 'legend.fontsize': 12})

    # Dictionaries for neater plot labels
    COSMO_DICT = {True: "Default cosmology ($h=0.673, \Omega_m=0.315$)", False: "Nondefault cosmology ($h=0.71, \Omega_m=0.27$)"}

    MAX_STEP_MPC = [50, 1, 1e-6]  # for 1.1e8 eV injection
    MAX_STEP_MPC = [50, 1, 0.01]  # for 1 TeV injection; 1e-6 step size is too computationally expensive for 1 Mpc trajectory
    E_INJ = 1e12
    MAX_TRAJ_MPC = [1]

    FN_E_INJ_DICT = {1.1e8 : "0.00011TeV",
                     1e10 : "0.001TeV",
                     1e12 : "1.0TeV", 
                    }

    for bp in [True]:
        for z in sim.REDSHIFTS:
            for i, max_traj in enumerate(MAX_TRAJ_MPC):
                # Want trajectories with the source at the same redshift and same max trajectory traced on the same plot
                fn_list_per_traj_and_z = []
                for max_step in MAX_STEP_MPC:
                    fn =  f"output/electrons/break_energy/bx/z{z}/inj1/delta/{FN_E_INJ_DICT[E_INJ]}/bamp1e-13G/{BP_DICT[bp]}/tol1e-09/iter0/set_cosmo/max_traj_{max_traj}Mpc/max_step_{max_step}Mpc/grid_5.0e+01Mpc/trajectory3d.txt"
                    fn_list_per_traj_and_z.append(fn)
                # Outside loop to plot stack of files for fixed trajectory
                print(fn_list_per_traj_and_z)
                subplot_binned_2d_trajectories(fn_data_list=fn_list_per_traj_and_z, 
                                               max_step_list=MAX_STEP_MPC, redshift=z, 
                                               traj=max_traj, 
                                               e_inj=E_INJ, 
                                               plot_prim=True,
                                               n=200)
                plot_energy_versus_x(fn_data_list=fn_list_per_traj_and_z,
                                     max_step_list=MAX_STEP_MPC, 
                                     redshift=z, 
                                     e_inj=E_INJ, 
                                     traj=max_traj, 
                                     use_bp=True)

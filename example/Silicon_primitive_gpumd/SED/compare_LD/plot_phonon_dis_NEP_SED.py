#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: LiangTing
2023/12/18 12:06:31
"""
from pylab import *
import numpy as np
import seaborn as sns

# for structure
from ase.io import read
from pySED.construct_BZ import BZPathHelper

#*************************** Set Seaborn style *************************
sns.set(style="ticks")
# Customize axis line, tick, and label properties
sns.set_context("paper", rc={"axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, 
    	                      "axes.labelsize": 18, "xtick.labelsize": 15.0, "ytick.labelsize": 15.0})
#*************************** Set Seaborn style *************************                

def get_path_and_distance():
	
    prim = read('POSCAR_prim')
    unitcell = prim.copy()
    dim = (20, 20, 20)
    supercell = unitcell.repeat(dim)

    lat = BZPathHelper(prim.cell, supercell.cell)
    q_path = np.array([[[0.0, 0.0, 0.0], [0.5, 0.0, 0.5]], [[0.5, 0.0, 0.5], [0.625, 0.25, 0.625]],
                       [[0.375, 0.375, 0.75], [0.0, 0.0, 0.0]], [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]]])

    file_path = "distances.txt"
    distances = read_distances_simple(file_path)
    
    pySED_paths, pySED_dists = [], []
    for p, d in zip(q_path, distances):
        p0, p1 = p[[0, -1]]

        pySED_path, pySED_dist = lat.build_commensurate_path(p0, p1)
        pySED_paths.append(pySED_path)
        d0, d1 = d[[0, -1]]
        pySED_dists.append(d0 + (d1 - d0) * pySED_dist)

    return pySED_paths, pySED_dists

def read_distances_simple(file_path):

    distances = []
    with open(file_path, 'r') as f:
        current_band = [] 
        for line in f:
            if line.startswith("# Band"):
                if current_band:
                    distances.append(np.array(current_band, dtype=float))
                    current_band = []  
            else:
                try:
                    current_band.append(float(line.strip()))
                except ValueError:
                    pass

        if current_band:
            distances.append(np.array(current_band, dtype=float))

    return distances         

def plot_dispersion_with_scaling(vmin=None, vmax=None):
    # Load NEP-LD data
    nep_data = np.loadtxt('band_structure.txt') 
    kpoint = np.loadtxt('qx_axis_data.txt')  
    
    # Load SED data
    out_files_name = 'silicon_300K'
    thz = np.loadtxt(out_files_name + '.THz')  
    sed_avg = np.loadtxt(out_files_name + '.SED')
    
    # Get SED distance 
    pySED_paths, pySED_dists = get_path_and_distance()
    #print(pySED_paths)

    # Plot
    if not vmin:
        vmin = np.trunc(sed_avg.min())
    if not vmax:
        vmax = np.trunc(sed_avg.max())

    fig, ax = plt.subplots(figsize=(6.5, 5))
    color = 'RdBu_r'
    sed_avg = np.log(sed_avg)
    levels = np.linspace(vmin, vmax, 350)

    # SED contour plot
    im = ax.contourf(np.hstack(pySED_dists), thz, sed_avg, cmap=color, levels=levels, vmin=vmin, vmax=vmax)

    # Colorbar
    ticks = np.arange(vmin, vmax + 0.1, 2)
    bar = fig.colorbar(im, ax=ax)
    bar.set_ticks(ticks)
    bar.set_ticklabels([f'{int(t)}' for t in ticks])
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, width=0, length=0, pad=0.6)
    bar.set_label(r'log($\Phi$($\mathbf{q}$, $\omega$)) (J $\cdot$ s)', fontsize=15)

    # Plot LD data
    lw = 1.5
    alpha = 0.8
    plot(nep_data[:, 0], nep_data[:, 1], color='grey', linestyle="--", lw=lw, alpha=alpha, label='NEP-LD')
    plot(nep_data[:, 0], nep_data[:, 2:], color='grey', linestyle="--", lw=lw, alpha=alpha)
    	
    for x in kpoint[1:-1]:  # Exclude first and last points
        ax.axvline(x, color='grey', linestyle='--', linewidth=0.8)

    # High symmetry points
    ax.set_xticks(kpoint)
    ax.set_xticklabels([r'$\Gamma$', 'X', 'U/K', r'$\Gamma$', 'L'], fontsize=17)

    # Axis settings
    xlim([0, kpoint[-1]])
    ylim([0, 16])
    ax.set_yticks(np.arange(0, 17, 4))
    ylabel('Frequency (THz)')

    savefig("Silicon.png", format='png', dpi=650, bbox_inches='tight')
    show()


if __name__ == "__main__":
    vmin = -25
    vmax = -9
    plot_dispersion_with_scaling(vmin, vmax)



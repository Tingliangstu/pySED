#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
@author: LiangTing
2023/12/18 12:06:31
"""
from pylab import *
import numpy as np
import seaborn as sns

#*************************** Set Seaborn style *************************
sns.set(style="ticks")
# Customize axis line, tick, and label properties
sns.set_context("paper", rc={"axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, 
    	                      "axes.labelsize": 18, "xtick.labelsize": 15.0, "ytick.labelsize": 15.0})
#*************************** Set Seaborn style *************************                

def plot_dispersion(vmin=None, vmax=None):

    ################################## load NEP-LD #############################
    nep_data = loadtxt('band_structure.txt')
    kpoint = loadtxt('qx_axis_data.txt')
    
    ################# load SED (depends on your own case) ######################
    out_files_name = 'graphene'
    thz = np.loadtxt(out_files_name + '.THz')
    sed_avg = np.loadtxt(out_files_name + '.SED')

    ## load SED-Q_distances_and_labels
    q_distances = []
    q_labels = {}
    with open(out_files_name + '.Q_distances_and_labels', 'r') as f:
        lines = f.readlines()
        # q_distances
        if lines[0].strip() == "Global distances along the paths:":
            q_distances = [float(x) for x in lines[1].split()]
        # q_labels
        if lines[3].strip() == "High-symmetry points and their distances:":
            for line in lines[4:]:
                distance, label = line.split(maxsplit=1)
                q_labels[float(distance)] = label.strip()
                
    # Scale q_distances to match kpoint
    scale_factor = kpoint[-1] / q_distances[-1]
    q_distances_scaled = [x * scale_factor for x in q_distances]
    
    # Scale q_labels keys
    q_labels_scaled = {k * scale_factor: v for k, v in q_labels.items()}
    
    ######################### Plot ##################
    # colarbar
    if not vmin:
        vmin = np.trunc(sed_avg.min())
    if not vmax:
        vmax = np.trunc(sed_avg.max())
    
    fig, ax = plt.subplots(figsize=(6.0, 5))
    
    color = 'RdBu_r'                        # 'jet', 'inferno'
    sed_avg = np.log(sed_avg)
    levels = np.linspace(vmin, vmax, 350)
    
    ax = gca()    
    im = ax.contourf(q_distances_scaled, thz, sed_avg, cmap=color, levels=levels, vmin=vmin, vmax=vmax)
    
        
    ticks = np.arange(vmin, vmax+0.01, 2)
    
    bar = fig.colorbar(im, ax=ax)
    bar.set_ticks(ticks)
    bar.set_ticklabels([f'{int(t)}' for t in ticks])
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, width=0, length=0, pad=0.6)
    bar.set_label(r'log($\Phi$($\mathbf{q}$, $\omega$)) (J $\cdot$ s)', fontsize=15)
    
    
    # For nep data
    lw = 1.5
    alpha = 0.6
    plot(nep_data[:, 0], nep_data[:, 1], color='grey', linestyle="--", lw=lw, alpha=alpha, label='NEP-LD')
    plot(nep_data[:, 0], nep_data[:, 2:], color='grey', linestyle="--", lw=lw, alpha=alpha)
    	
    for x in kpoint[1:-1]:  # Exclude first and last points
        ax.axvline(x, color='grey', linestyle='--', linewidth=0.8)

    # labels set
    #xticks = list(q_labels_scaled.keys())
    #xticklabels = [r'$\Gamma$' if label == 'G' else label for label in q_labels_scaled.values()]  # Replace 'G' with '$\Gamma$'
    
    gca().set_xticks(kpoint)
    gca().set_xticklabels([r'$\Gamma$', 'M', 'K', r'$\Gamma$'], fontsize=16)
    
    xlim([0, max(kpoint)])
    ylim([0, 50])
    
    ylabel('Frequency (THz)')
    legend(loc="best", frameon=False, fontsize=15)

    savefig("Graphene.png", format='png', dpi=650, bbox_inches='tight')
    show()


if __name__ == "__main__":
	
	  # bar min and max
    vmin = -26
    vmax = -8
    plot_dispersion(vmin, vmax)

    print('******************** All Done !!! *************************')

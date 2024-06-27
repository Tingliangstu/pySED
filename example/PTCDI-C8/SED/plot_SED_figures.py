# -*- coding: utf-8 -*-
"""
Created on 2022/10/10 15:43:40
@author: Liang Ting
"""
import numpy as np
import matplotlib.pyplot as plt

def plot_sed(out_files_name, plot_cutoff_freq = None, plot_interval = 5, vmin = None, vmax = None, if_show_figures=False):
    # load data for draw
    sed_avg = np.loadtxt(out_files_name + '.SED')
    qpoints = np.loadtxt(out_files_name + '.Qpts')
    thz = np.loadtxt(out_files_name + '.THz')

    # ******************** Control plotting params ********************
    log = True
    color = 'inferno'       # 'inferno', 'jet', 'Spectral', 'RdYlGn' (https://matplotlib.org/stable/tutorials/colors/colormaps.html)

    interp = 'hanning'       # 'hanning'
    df = plot_interval       # Scale interval for drawing
    fontsize = 20

    # For plot scale
    max_thz = max(thz)
    max_sed_y = np.size(sed_avg, 0)
    scale_factor = max_sed_y / max_thz

    # ******************** Whether to apply log scaling data ********************
    if log:
        sed_avg = np.log(sed_avg)

    ### ******************** Creat a figure, set its size ********************
    fig, ax = plt.subplots()
    fig.set_size_inches(7, 9, forward=True)       # Control the size of the output image
    fig.tight_layout(pad=6)
    if not vmin:
        vmin = np.trunc(sed_avg.min())
    if not vmax:
        vmax = np.trunc(sed_avg.max())

    print('********* The vmin = {} and vmax = {} *********'.format(vmin, vmax))
    ticks = np.arange(vmin, vmax, 1)

    # Plot colorbar
    plt.rcParams['font.size'] = 16
    im = ax.imshow(sed_avg, cmap = color, interpolation = interp, aspect = 'auto', origin = 'lower', vmax=vmax, vmin=vmin)
    bar = fig.colorbar(im, ax=ax, label = r'log($\Phi $($\bfq$, $\omega)$) (arbitrary unit)')
    bar.set_ticks(ticks)
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=10, width=0, length=0, pad=5)

    # configure the plot axix line
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2.5)

    # yticks
    freqs = np.arange(0, thz.max(), df)
    nf = len(freqs)
    ids = np.zeros(nf)
    for i in range(nf):
        ids[i] = np.argwhere(thz <= freqs[i]).max()
    ax.set_yticks(ids)
    ax.set_yticklabels(list(map(str, freqs)), fontsize=16)
    ax.tick_params(which='both', axis='y', width=2)
    ax.tick_params(which='major', length=5)

    # xticks
    xticks = [-0.5, len(qpoints)-0.5]
    ax.set_xticks(xticks)
    ax.set_xticklabels([r'$\Gamma$', r'X'], fontsize=18)

    # ******************** set the figure labels ********************
    ax.set_xlabel(r'$\bfq$',labelpad = 5.0,fontweight='normal',fontsize=fontsize)
    ax.set_ylabel('Frequency (THz)', labelpad = 8.0, fontweight='normal', fontsize=fontsize)
    ax.tick_params(bottom=False, top=False, left=True, right=False)

    fig.suptitle('PTCDI-C8 (300 K)', x=0.48, y=0.95, fontweight='normal', fontsize=fontsize)

    if plot_cutoff_freq:
        ax.set_ylim([0, plot_cutoff_freq * scale_factor])

    plt.savefig('{}-SED-nice.png'.format(out_files_name), format='png', dpi = 600, bbox_inches='tight')

    if if_show_figures:
        plt.show()

if __name__ == "__main__":

    out_files_name = 'PTCDI-C8'
    plot_cutoff_freq = 1               # THz
    plot_interval = 0.5                # THz
    vmax = -8                          # For SED image (high)
    vmin = -20                         # For SED image (low)

    plot_sed(out_files_name, plot_cutoff_freq, plot_interval=plot_interval, vmax=vmax, vmin=vmin)

    print('************** ALL DONE !!! **************')

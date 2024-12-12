'''
@author:
**************************  LiangTing ****************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2024/12/6 23:03:21 **********************
'''
from pylab import *
import seaborn as sns
import numpy as np

# *************************** Set Seaborn style for good look*************************
sns.set(style="ticks")
# Customize axis line, tick, and label properties
sns.set_context("paper", rc={"axes.linewidth": 0.75, "xtick.major.width": 0.75, "ytick.major.width": 0.75,
                             "axes.labelsize": 16, "xtick.labelsize": 13.5, "ytick.labelsize": 13.5})
# *************************** Set Seaborn style *************************

def plot_bands(data, params):
    # Get data
    sed_avg = data.sed_avg
    qpoints = data.qpoints
    thz = data.freq_fft
    q_distances = data.q_distances
    q_labels = data.q_labels

    # ******************** Control plotting params ********************
    color = 'RdBu_r'                        # 'jet', 'inferno'
    interp = 'hanning'                      # 'hanning'
    df = params.plot_interval               # Scale interval for drawing

    # For plot scale
    max_thz = max(thz)
    max_sed_y = np.size(sed_avg, 0)
    scale_factor = max_sed_y / max_thz
    vmin = None   # -8 for test
    vmax = None   # -20 for test

    # ******************** Whether to apply log scaling data ********************
    sed_avg = np.log(sed_avg)
    ### ******************** Creat a figure, set its size ********************

    fig, ax = plt.subplots()
    fig.set_size_inches(5.5+(params.num_qpaths/2-0.5), 5)  # Control the size of the output image

    if not vmin:
        vmin = np.trunc(sed_avg.min())
    if not vmax:
        vmax = np.trunc(sed_avg.max())

    levels = np.linspace(vmin, vmax, 800)

    """
    #Choose the plotting method based on the number of qpaths
    Set yticks (This part of the code is redundant, maybe I will remove the params.num_qpaths == 1,
     for now just to better repeat the previous results.)
    """
    if params.num_qpaths > 1 or params.use_contourf:
        im = ax.contourf(q_distances, thz, sed_avg, cmap=color, levels=levels, vmin=vmin, vmax=vmax)
        # Draw vertical grey lines for q_labels (excluding endpoints)
        keys = list(q_labels.keys())
        for x in keys[1:-1]:            # Exclude first and last points
            ax.axvline(x, color='grey', linestyle='--', linewidth=0.8)

    else:
        print("You are using \'imshow\' method for SED plotting since single Qpaths, "
              "it can be change to \'contourf\' by setting \'use_contourf = 1\' in the input.in.")
        im = ax.imshow(sed_avg, cmap=color, interpolation=interp, aspect='auto', origin='lower', vmax=vmax, vmin=vmin)

    # colarbar
    ticks = np.arange(vmin, vmax+0.01, 2)
    bar = fig.colorbar(im, ax=ax)
    bar.set_ticks(ticks)
    bar.set_ticklabels([str(int(t)) for t in ticks])
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, width=0, length=0, pad=0.6)
    bar.set_label(r'log($\Phi$($\mathbf{q}$, $\omega$)) (J $\cdot$ s)', fontsize=13.5)

    # Set xticks and xticklabels
    if params.num_qpaths > 1 or params.use_contourf:
        xticks = list(q_labels.keys())

    else:
        ax.set_xlim([0, len(qpoints)-1])
        xticks = [0, len(qpoints)-1]

    xticklabels = [r'$\Gamma$' if label == 'G' else label for label in q_labels.values()]  # Replace 'G' with '$\Gamma$'
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=16)

    """
    Set yticks (This part of the code is redundant, maybe I will remove the params.num_qpaths == 1,
     for now just to better repeat the previous results.)
    """
    if params.num_qpaths > 1 or params.use_contourf:
        max_freq = params.plot_cutoff_freq if params.plot_cutoff_freq else thz.max()
        num_ticks = int(np.ceil(max_freq / df)) + 1
        freqs = np.linspace(0, max_freq, num_ticks)

        ax.set_yticks(freqs)
        yticks_labels = [f'{f:.1f}'.rstrip('0').rstrip('.') for f in freqs]
        ax.set_yticklabels(yticks_labels, fontsize=13.5)
        ax.set_ylim([0, max_freq])

    else:
        freqs = np.arange(0, np.ceil(thz.max()) + 0.01, df)
        ids = np.zeros(len(freqs))
        for i in range(len(ids)):
            ids[i] = np.argwhere(thz <= freqs[i]).max()

        ax.set_yticks(ids)
        yticks_labels = [f'{f:.1f}'.rstrip('0').rstrip('.') for f in freqs]
        ax.set_yticklabels(yticks_labels, fontsize=13.5)

        # ax.set_ylim()
        if params.plot_cutoff_freq:
            ax.set_ylim([0, params.plot_cutoff_freq * scale_factor])

    ax.set_ylabel('Frequency (THz)', fontsize=16)

    plt.savefig('{}-SED.png'.format(params.out_files_name), format='png', dpi=650, bbox_inches='tight')

    if params.if_show_figures:
        plt.show()

def lorentzian(xarr, center, amplitude, hwhm):
    return amplitude / (1 + ((xarr - center) / hwhm) ** 2)

def plot_slice(data, params):
    # ******************** Get data ********************
    sed_avg = data.sed_avg
    qpoints = data.qpoints
    thz = data.freq_fft
    q_index = params.q_slice_index

    ### ******************** creat a figure, set its size ********************
    fig, ax = plt.subplots()
    fig.set_size_inches(8, 4)

    # ******************** For saving different files ********************
    save_flag = False
    # Plot
    alpha = 0.6

    ax.semilogy(thz, sed_avg[:, q_index], ls='-', lw=1.4, color="#aa3474", marker='o', ms=5.5, fillstyle='full', alpha=alpha)

    if params.plot_lorentz:
        save_flag = True
        total = np.zeros(len(sed_avg[:, q_index]))
        for i in range(len(params.popt[:, 0])):
            if params.popt[i, 2] == 0:
                continue
            ax.semilogy(thz, lorentzian(thz, params.popt[i, 0], params.popt[i, 1], params.popt[i, 2]),
                        ls='-', lw=1.5, color='C2', alpha=alpha)

            total = total + (lorentzian(thz, params.popt[i, 0], params.popt[i, 1], params.popt[i, 2]))

        ax.semilogy(thz, total, ls='--', lw=1.5, color='grey', alpha=alpha+0.2)

    # ******************** set the figure labels ********************
    ax.set_ylabel(r'log($\Phi$($\omega)$) (J $\cdot$ s)')
    ax.set_xlabel('Frequency (THz)')
    fig.suptitle(r'$\mathbf{{q}}$ = ({0:.3f}, {1:.3f}, {2:.3f})'.format(qpoints[q_index, 0], qpoints[q_index, 1],
                                                               qpoints[q_index, 2]), y=0.95, fontsize=15)

    if params.lorentz_fit_cutoff:
        ax.set_xlim([0, params.lorentz_fit_cutoff])
    else:
        ax.set_xlim([0, np.ceil(thz.max())+0.01])

    if save_flag:
        plt.savefig('LORENTZ-fitting-{}-qpoint.png'.format(params.q_slice_index), format='png', dpi=650, bbox_inches='tight')

    else:
        plt.savefig('SED-{}-qpoint.png'.format(params.q_slice_index), format='png', dpi=650, bbox_inches='tight')

    if params.if_show_figures:
        plt.show()

    plt.close(fig)

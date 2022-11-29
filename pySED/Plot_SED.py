'''
@author:
**************************  LiangTing ****************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/5/15 23:03:21 **********************
'''

import numpy as np
import matplotlib.pyplot as plt

def plot_bands(data, params):
    # Get data
    sed_avg = data.sed_avg
    qpoints = data.qpoints
    thz = data.freq_fft

    # ******************** Control plotting params ********************
    log = True
    color = 'jet'                           # 'inferno'
    interp = 'hanning'                      # 'hanning'
    df = params.plot_interval               # Scale interval for drawing

    # For plot scale
    max_thz = max(thz)
    max_sed_y = np.size(sed_avg, 0)
    scale_factor = max_sed_y / max_thz
     
    # ******************** Whether to apply log scaling data ********************
    if log:
        sed_avg = np.log(sed_avg)

    ### ******************** Creat a figure, set its size ********************
    fig, ax = plt.subplots()
    fig.set_size_inches(4, 6, forward=True)  # Control the size of the output image
    fig.tight_layout(pad=5)                  #
    vmin = np.trunc(sed_avg.min())
    vmax = np.trunc(sed_avg.max())

    im = ax.imshow(sed_avg, cmap = color, interpolation = interp, aspect = 'auto', origin = 'lower', vmax=vmax, vmin=vmin)

    fig.colorbar(im, ax=ax, label = r'log($\Phi$($\omega)$)  (J.s)')
    # configure the plot
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

    freqs = np.arange(0, thz.max(), df)
    nf = len(freqs)
    ids = np.zeros(nf)
    for i in range(nf):
        ids[i] = np.argwhere(thz <= freqs[i]).max()
    ax.set_yticks(ids)
    ax.set_yticklabels(list(map(str, freqs)))

    xticks = [0, len(qpoints)-1]
    ax.set_xticks(xticks)
    xlabels = ['']*len(xticks)
    for i in range(len(xticks)):
        xlabels[i] = '({:.1f},{:.1f},{:.1f})'.format(qpoints[xticks[i], 0],
                qpoints[xticks[i], 1], qpoints[xticks[i], 2])

    ax.set_xticklabels(xlabels)

    ax.minorticks_on()
    ax.tick_params(which = 'both', axis='y', width=1, labelsize='large')
    ax.tick_params(which = 'both', axis='x', width=1, labelsize='large',
            labelrotation = 0.0,pad = 5.0)
    ax.tick_params(which = 'major', length = 5)
    ax.tick_params(which = 'minor', length = 3, color='k')
    # plt.tick_params(axis = 'x', which = 'both', labelbottom = False)

    # ******************** set the figure labels ********************
    ax.set_xlabel(r'$\bfq$',labelpad = 5.0,fontweight='normal',fontsize='x-large')
    ax.set_ylabel(r'Frequency (THz)',labelpad = 3.0,fontweight='normal',fontsize='x-large')

    fig.suptitle(r'$\Phi$($\bfq$,$\omega)$',y = 0.95,fontsize='x-large')

    # ax.set_xlim()
    if params.plot_cutoff_freq:
        ax.set_ylim([0, params.plot_cutoff_freq * scale_factor])

    plt.savefig('{}-SED.png'.format(params.out_files_name), format='png', dpi = 600, bbox_inches='tight')

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

    # ******************** Whether to apply log scaling data ********************
    log = True

    df = 5                                              # Scale interval for drawing
    freqs = np.arange(0, thz.max(), df)
    nf = len(freqs)
    ids = np.zeros(nf)

    for i in range(nf):
        ids[i] = np.argwhere(thz <= freqs[i]).max()

    ### ******************** creat a figure, set its size ********************
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 6, forward=True)
    fig.tight_layout(pad=10)

    # ******************** Configure the plot ********************
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)

    # ******************** For saving different files ********************
    save_flag = False

    if log:
        ax.semilogy(thz, sed_avg[:, q_index], ls='-', lw=1.5, color='k', marker = 'x', ms =3, mec='lightcoral')
        if params.plot_lorentz:
            save_flag = True
            total = np.zeros(len(sed_avg[:, q_index]))
            for i in range(len(params.popt[:, 0])):
                if params.popt[i, 2] == 0:
                    continue
                ax.semilogy(thz, lorentzian(thz,params.popt[i, 0], params.popt[i, 1], params.popt[i, 2]),
                    ls = '-', lw = 2, color = 'm')

                total = total + (lorentzian(thz, params.popt[i, 0], params.popt[i, 1], params.popt[i, 2]))

            ax.semilogy(thz, total, ls = '--', lw = 1, color = 'k')

    else:
        ax.plot(sed_avg[:, q_index], ls='-', lw=1, color='k', marker='o', ms=2, mfc='b', mec='k', mew=1)

        if params.plot_lorentz:
            total = np.zeros(len(sed_avg[:, q_index]))
            for i in range(len(params.popt[:, 0])):
                if params.popt[i, 2] == 0:
                    continue
                ax.plot(lorentzian(np.arange(len(sed_avg[:, q_index])),
                                   params.popt[i, 0], params.popt[i, 1], params.popt[i, 2]),
                        ls='-', lw=1, marker='o', mfc='r', mec='r', ms=1, mew=0, color='r')
                total = total + (lorentzian(np.arange(len(sed_avg[:, q_index])),
                                            params.popt[i, 0], params.popt[i, 1], params.popt[i, 2]))
            ax.plot(total, ls='-', lw=0.5, marker='o', mfc='b', mec='b',
                    ms=0.5, mew=0, color='b')

    ax.minorticks_on()
    ax.tick_params(which='both', width=1, labelsize = 'x-large')
    ax.tick_params(which='major', length=5)
    ax.tick_params(which='minor', length=3, color='k')
    # plt.tick_params(axis='x',which='both',labelbottom=False)

    # ******************** set the figure labels ********************
    ax.set_ylabel(r'log($\Phi$($\omega)$)', labelpad=35.0, fontweight='normal',
                  fontsize='x-large', rotation='horizontal')
    ax.set_xlabel('Frequency (THz)', labelpad=3.0, fontweight='normal', fontsize='x-large')

    # ax.set_xlabel(r'$\omega$ (THz)', labelpad=3.0, fontweight='normal', fontsize='large')

    fig.suptitle(r'$\bfq$ = ({:.3f}, {:.3f}, {:.3f})'.format(
        qpoints[q_index, 0], qpoints[q_index, 1], qpoints[q_index, 2]), y=0.80, fontsize='x-large')

    if params.lorentz_fit_cutoff:
        ax.set_xlim([0, params.lorentz_fit_cutoff])

    if save_flag:
        plt.savefig('LORENTZ-fitting-{}-qpoint.png'.format(params.q_slice_index), format='png',dpi = 600, bbox_inches='tight')
    else:
        plt.savefig('SED-{}-qpoint.png'.format(params.q_slice_index), format='png', dpi=600, bbox_inches='tight')

    if params.if_show_figures:
        plt.show()

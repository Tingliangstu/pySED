# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED.  If not, see <http://www.gnu.org/licenses/>.
# =============================================================================

"""Plot saved SED maps, single-q spectra, fit overlays, and lifetimes.

The data object supplies spectra and axes from FileIO.load_data. Its sed_unit
controls intensity labels; changing an input flag does not reinterpret a
marked file. Plotting uses logarithmic intensity display without rescaling
the numerical spectrum used by the fitter.
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

from pySED.postprocess import output_files

# Use one Seaborn theme for every pySED figure in this plotting module.
sns.set_theme(
    context="paper",
    style="ticks",
    rc={
        "axes.grid": False,
        "axes.linewidth": 0.75,
        "axes.labelsize": 16,
        "xtick.major.width": 0.75,
        "ytick.major.width": 0.75,
        "xtick.labelsize": 13.5,
        "ytick.labelsize": 13.5,
    },
)

def _resolve_color_limits(sed, params):
    """Return natural-log color limits in the spectrum's stored output scale."""

    values = np.asarray(sed)
    positive = values[np.isfinite(values) & (values > 0)]
    if positive.size == 0:
        raise ValueError("SED data contain no positive values to plot")

    log_values = np.log(positive)
    automatic_min = np.trunc(log_values.min())
    automatic_max = np.trunc(log_values.max())
    if np.isclose(automatic_min, automatic_max):
        automatic_min -= 1.0
        automatic_max += 1.0

    vmin = (
        params.colorbar_min
        if params.colorbar_min is not None
        else automatic_min
    )
    vmax = (
        params.colorbar_max
        if params.colorbar_max is not None
        else automatic_max
    )
    if vmin >= vmax:
        raise ValueError(
            "resolved colorbar_min must be smaller than colorbar_max"
        )
    return vmin, vmax


def _prepare_sed_for_log_scale(sed, vmin, vmax):
    """Take the natural log of positive SED values and clip for rendering."""

    values = np.asarray(sed)
    positive = np.ma.masked_less_equal(values, 0.0)
    finite = np.ma.masked_invalid(positive)
    return np.ma.clip(np.ma.log(finite), vmin, vmax)


def _print_plot_detail(label, value):
    """Print one plot detail using the pipeline's column alignment."""

    print(f"      {label:<16}: {value}")


def _input_filename(params):
    """Return the actual input filename for user-facing plot guidance."""

    return os.path.basename(getattr(params, "input_file", "input_SED.in"))


def _print_color_scale(params, vmin, vmax, unit="eV/THz"):
    """Print the resolved color scale in aligned, readable rows."""

    min_source = "user-defined" if params.colorbar_min is not None else "automatic"
    max_source = "user-defined" if params.colorbar_max is not None else "automatic"
    input_name = _input_filename(params)
    _print_plot_detail(
        "Colorbar",
        f"ln[SED / ({unit})] (dimensionless)",
    )
    _print_plot_detail(
        "colorbar_min",
        f"{vmin:.6g} ({min_source}, ln scale)",
    )
    _print_plot_detail(
        "colorbar_max",
        f"{vmax:.6g} ({max_source}, ln scale)",
    )
    _print_plot_detail(
        "Fine-tune with",
        f"colorbar_min / colorbar_max in {input_name}",
    )


def resolve_slice_frequency_limit(thz, params, *, fitting=False):
    """Return the displayed upper frequency and its input control."""

    lorentz_freq_max = getattr(params, 'lorentz_fit_freq_max', None)
    if (
        fitting
        and lorentz_freq_max is not None
    ):
        return lorentz_freq_max, "lorentz_fit_freq_max"
    plot_cutoff = getattr(params, 'plot_cutoff_freq', None)
    if plot_cutoff is not None:
        return plot_cutoff, "plot_cutoff_freq"
    return float(np.max(thz)), None


def resolve_slice_frequency_start(thz, params, *, fitting=False):
    """Return the displayed lower frequency and its input control."""

    lorentz_freq_min = getattr(params, 'lorentz_fit_freq_min', None)
    if (
        fitting
        and lorentz_freq_min is not None
    ):
        return lorentz_freq_min, "lorentz_fit_freq_min"
    return float(np.min(thz)), None


def plot_bands(data, params):
    """Write a dispersion heatmap from SED values shaped (frequency, q).

    freq_fft supplies THz frequencies; q_distances and q_labels describe the
    path. Multiple segments or use_contourf select physical-coordinate contours;
    otherwise imshow displays the matrix by pixel index. Color limits apply to
    the natural logarithm of the stored intensities. The colorbar names the
    scale recorded in data, and output_files selects the output filename.
    """
    sed_avg = data.sed_avg
    qpoints = data.qpoints
    thz = data.freq_fft
    q_distances = data.q_distances
    q_labels = data.q_labels

    color = params.plot_color               # 'inferno', 'RdBu_r', 'jet', 'Spectral', 'RdYlGn' (https://matplotlib.org/stable/tutorials/colors/colormaps.html)
    interp = 'hanning'                      # 'hanning'
    df = params.plot_interval               # Scale interval for drawing

    max_thz = max(thz)
    max_sed_y = np.size(sed_avg, 0)
    scale_factor = max_sed_y / max_thz

    vmin, vmax = _resolve_color_limits(sed_avg, params)
    sed_for_plot = _prepare_sed_for_log_scale(sed_avg, vmin, vmax)


    fig, ax = plt.subplots()
    fig.set_size_inches(5.5+(params.num_qpaths/2-0.5), 5)  # Control the size of the output image

    unit = getattr(data, 'sed_unit', 'eV/THz')
    _print_color_scale(params, vmin, vmax, unit)

    # Contour coordinates follow physical Q distances and frequencies;
    # imshow uses pixel indices. Keep their coordinate handling distinct.
    if params.num_qpaths > 1 or params.use_contourf:
        levels = np.linspace(vmin, vmax, 350)

        im = ax.contourf(
            q_distances,
            thz,
            sed_for_plot,
            cmap=color,
            levels=levels,
            vmin=vmin,
            vmax=vmax,
        )

        # Draw vertical grey lines for q_labels (excluding endpoints)
        keys = list(q_labels.keys())
        for x in keys[1:-1]:            # Exclude first and last points
            ax.axvline(x, color='grey', linestyle='--', linewidth=0.8)

    else:
        _print_plot_detail(
            "Rendering",
            "SED is rendered with imshow; set "
            f"use_contourf = 1 in {_input_filename(params)} "
            "to use contourf",
        )
        im = ax.imshow(
            sed_for_plot,
            cmap=color,
            interpolation=interp,
            aspect='auto',
            origin='lower',
            vmin=vmin,
            vmax=vmax,
        )

    bar = fig.colorbar(im, ax=ax)
    ticks = np.arange(vmin, vmax + 0.01, 2.0)
    if ticks.size >= 2:
        bar.set_ticks(ticks)
        bar.set_ticklabels([f"{tick:g}" for tick in ticks])
    bar.outline.set_visible(False)
    bar.ax.tick_params(labelsize=8, width=0, length=0, pad=0.6)
    colorbar_label = (
        r'$\ln\!\left[\mathrm{SED}(\mathbf{q},\omega)'
        r'/(\mathrm{eV/THz})\right]$'
        if unit == 'eV/THz' else f'ln[SED / ({unit})]'
    )
    bar.set_label(colorbar_label, fontsize=13.5)

    # Set xticks and xticklabels
    if params.num_qpaths > 1 or params.use_contourf:
        xticks = list(q_labels.keys())

    else:
        ax.set_xlim([0, len(qpoints)-1])
        xticks = [0, len(qpoints)-1]

    xticklabels = [r'$\Gamma$' if label == 'G' else label for label in q_labels.values()]  # Replace 'G' with '$\Gamma$'
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=16)

    # Map physical frequency labels to the selected renderer's coordinates.
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

        if params.plot_cutoff_freq:
            ax.set_ylim([0, params.plot_cutoff_freq * scale_factor])

    ax.set_ylabel('Frequency (THz)', fontsize=16)

    output_path = output_files.dispersion_figure(params)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(output_path), format='png', dpi=650, bbox_inches='tight')

    if params.if_show_figures:
        plt.show()

def plot_slice(data, params, fit_result=None):
    """Write one q-point intensity spectrum, optionally with fitted curves.

    params.qpoint_slice_index selects a column of data.sed_avg. Plot original
    intensities on a logarithmic y axis against freq_fft in THz. When a
    fit_result is supplied, overlay its peak_fits curves and use the fitting
    frequency bounds and filename; otherwise use the ordinary slice controls.
    The y label follows data.sed_unit. The figure is closed after saving.
    """
    sed_avg = data.sed_avg
    qpoints = data.qpoints
    thz = data.freq_fft
    q_index = params.qpoint_slice_index

    fig, ax = plt.subplots()
    fig.set_size_inches(8, 4)

    alpha = 0.6

    ax.semilogy(
        thz,
        sed_avg[:, q_index],
        ls='-',
        lw=1.4,
        color="#aa3474",
        marker='o',
        ms=5.5,
        fillstyle='full',
        alpha=alpha,
        label='SED data',
        zorder=5,
    )

    if fit_result is not None:
        peak_fits = fit_result.peak_fits
        fit_styles = {
            'lorentz': ('#45a65c', 'Lorentz peak fits'),
            'dho': ('#7f7f7f', 'DHO peak fits'),
        }
        shown_models = set()
        for fit in peak_fits:
            fit_model = fit.get('fitting_function', 'lorentz')
            fit_color, fit_label = fit_styles[fit_model]
            fit_frequency = np.asarray(fit['frequency'])
            fitted_curve = np.asarray(fit['predicted'])
            physical = fit_frequency >= 0.0
            ax.semilogy(
                fit_frequency[physical],
                fitted_curve[physical],
                ls='-',
                lw=2.2,
                color=fit_color,
                alpha=0.8,
                label=(
                    fit_label if fit_model not in shown_models else None
                ),
                zorder=7,
            )
            shown_models.add(fit_model)
        if peak_fits:
            handles, labels = ax.get_legend_handles_labels()
            by_label = dict(zip(labels, handles))
            legend_labels = ['SED data']
            for model in ('lorentz', 'dho'):
                label = fit_styles[model][1]
                if label in by_label:
                    legend_labels.append(label)
            ax.legend(
                [by_label[label] for label in legend_labels],
                legend_labels,
                fontsize=13,
                frameon=False,
                handlelength=2.8,
                labelspacing=0.6,
            )

    ax.set_ylabel(
        r'$\Phi(\mathbf{q},\omega)$' + f" ({getattr(data, 'sed_unit', 'eV/THz')})",
        fontsize=16,
    )
    ax.set_xlabel('Frequency (THz)', fontsize=16)
    ax.tick_params(axis='both', labelsize=14)
    fig.suptitle(r'$\mathbf{{q}}$ = ({0:.3f}, {1:.3f}, {2:.3f})'.format(qpoints[q_index, 0], qpoints[q_index, 1],
                                                               qpoints[q_index, 2]), y=0.95, fontsize=15)

    min_frequency, _ = resolve_slice_frequency_start(thz, params, fitting=fit_result is not None)
    max_frequency, _ = resolve_slice_frequency_limit(thz, params, fitting=fit_result is not None)
    ax.set_xlim([min_frequency, max_frequency])

    output_path = str(output_files.slice_figure(params, fitting=fit_result is not None))
    output_directory = os.path.dirname(output_path)
    if output_directory:
        os.makedirs(output_directory, exist_ok=True)
    plt.savefig(
        output_path,
        format='png',
        dpi=650,
        bbox_inches='tight',
    )

    if params.if_show_figures:
        plt.show()

    plt.close(fig)


def plot_lifetime_summary(frequency, lifetime, params=None):
    """Write a frequency-versus-lifetime scatter plot and return its path.

    Inputs contain frequencies in THz and lifetimes in ps. Filter non-finite
    values, negative frequencies, and non-positive lifetimes; use a logarithmic
    lifetime axis. Raise ValueError if no valid points remain. The output goes
    beside total or partial SED figures according to params.
    """

    frequency = np.asarray(frequency, dtype=float)
    lifetime = np.asarray(lifetime, dtype=float)
    valid = (
        np.isfinite(frequency)
        & np.isfinite(lifetime)
        & (frequency >= 0.0)
        & (lifetime > 0.0)
    )
    if not np.any(valid):
        raise ValueError("no positive finite lifetimes are available to plot")

    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    sns.scatterplot(
        x=frequency[valid],
        y=lifetime[valid],
        ax=ax,
        s=30,
        color="#aa3474",
        alpha=0.72,
        linewidth=0,
        zorder=3,
    )
    ax.set_xlabel("Frequency (THz)", fontsize=16)
    ax.set_ylabel("Lifetime", fontsize=16)
    ax.set_yscale("log")
    ax.tick_params(axis="both", labelsize=14)
    ax.grid(False, which="both", axis="both")

    output_path = output_files.lifetime_summary_figure(params)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        format="png",
        dpi=450,
        bbox_inches="tight",
    )
    plt.close(fig)
    return str(output_path)

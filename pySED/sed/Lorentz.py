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

"""Independent zero-background fitting of SED spectral peaks.

The fitter deliberately keeps the pySED workflow simple. Peaks are detected
with a local-noise significance threshold, and each peak receives the complete
local basin between the smoothed valleys separating it from neighboring
detected peaks. Every peak is then fitted independently. The only optional
modelling choice is the line shape: Lorentzian, velocity-spectrum DHO, or an
AICc-based choice between those two shapes.
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, windows

from pySED.sed import FileIO, Plot_SED


_PEAK_DETECTION_SMOOTH_WINDOW = 7
_PEAK_DETECTION_NOISE_WINDOW = 31
_MAD_TO_SIGMA = 1.4826
_FITTING_FUNCTIONS = ("lorentz", "dho")
_AUTO_AICC_EQUIVALENCE = 2.0
_DHO_MIN_DISTINGUISHABLE_DAMPING_RATIO = 0.05


def _odd_window_length(size, requested):
    """Return the largest usable odd window, or one for very short data."""

    window_len = min(size, requested)
    if window_len % 2 == 0:
        window_len -= 1
    return max(window_len, 1)


def _hann_smooth(values, window_len=_PEAK_DETECTION_SMOOTH_WINDOW):
    """Return a detection-only Hann-smoothed copy of one spectrum."""

    values = np.asarray(values, dtype=float)
    window_len = _odd_window_length(values.size, window_len)
    if window_len < 3:
        return values.copy()

    window = windows.hann(window_len)
    pad = window_len // 2
    padded = np.pad(values, (pad, pad), mode="reflect")
    return np.convolve(padded, window / window.sum(), mode="valid")


def _rolling_mad(values, window_len=_PEAK_DETECTION_NOISE_WINDOW):
    """Return a rolling median absolute deviation with reflected edges."""

    values = np.asarray(values, dtype=float)
    window_len = _odd_window_length(values.size, window_len)
    if window_len < 3:
        return np.zeros_like(values)

    pad = window_len // 2
    padded = np.pad(values, (pad, pad), mode="reflect")
    local_windows = np.lib.stride_tricks.sliding_window_view(
        padded,
        window_len,
    )
    local_medians = np.median(local_windows, axis=1)
    return np.median(
        np.abs(local_windows - local_medians[:, np.newaxis]),
        axis=1,
    )


def _smoothed_log_sed(sed):
    """Return the seven-point detection-smoothed logarithmic spectrum."""

    sed = np.asarray(sed, dtype=float)
    positive = sed[np.isfinite(sed) & (sed > 0)]
    if not positive.size:
        raise ValueError(
            "peak detection requires positive finite SED values"
        )
    log_floor = max(np.min(positive) * 1.0e-12, np.finfo(float).tiny)
    log_sed = np.log(
        np.where(np.isfinite(sed) & (sed > 0), sed, log_floor)
    )
    return _hann_smooth(log_sed)


def _find_sed_peaks(sed, peak_min_significance=4.0):
    """Detect SED peaks using dimensionless local-noise significance."""

    sed = np.asarray(sed, dtype=float)
    if (
        not np.isfinite(peak_min_significance)
        or peak_min_significance <= 0
    ):
        raise ValueError("peak_min_significance must be positive and finite")

    smooth_log_sed = _smoothed_log_sed(sed)
    positive = sed[np.isfinite(sed) & (sed > 0)]
    log_floor = max(np.min(positive) * 1.0e-12, np.finfo(float).tiny)
    log_sed = np.log(
        np.where(np.isfinite(sed) & (sed > 0), sed, log_floor)
    )
    residual = log_sed - smooth_log_sed

    local_sigma = _MAD_TO_SIGMA * _rolling_mad(residual)
    residual_median = np.median(residual)
    global_sigma = _MAD_TO_SIGMA * np.median(
        np.abs(residual - residual_median)
    )
    sigma_floor = max(0.1 * global_sigma, np.finfo(float).eps)
    local_sigma = np.maximum(local_sigma, sigma_floor)

    minimum_prominence = peak_min_significance * local_sigma
    candidates, detection_properties = find_peaks(
        smooth_log_sed,
        prominence=minimum_prominence,
    )
    detection_significance = (
        detection_properties["prominences"] / local_sigma[candidates]
    )

    refinement_radius = _odd_window_length(
        sed.size,
        _PEAK_DETECTION_SMOOTH_WINDOW,
    ) // 2
    refined = {}
    for candidate, significance in zip(candidates, detection_significance):
        left = max(0, candidate - refinement_radius)
        right = min(sed.size, candidate + refinement_radius + 1)
        peak = left + int(np.argmax(sed[left:right]))
        if 0 < peak < sed.size - 1:
            refined[peak] = max(refined.get(peak, 0.0), significance)

    peaks = np.array(sorted(refined), dtype=int)
    significance = np.array([refined[peak] for peak in peaks], dtype=float)
    return peaks, significance


def _select_lorentz_frequency_range(
    frequency,
    sed,
    freq_min=None,
    freq_max=None,
):
    """Select the inclusive frequency interval used for fitting."""

    frequency = np.asarray(frequency)
    sed = np.asarray(sed)
    mask = np.ones(frequency.shape, dtype=bool)
    if freq_min is not None:
        mask &= frequency >= freq_min
    if freq_max is not None:
        mask &= frequency <= freq_max
    if np.count_nonzero(mask) < 3:
        raise ValueError(
            "the requested spectral fitting frequency range contains "
            "fewer than three frequency samples"
        )
    return frequency[mask], sed[mask]


def _lorentzian(xarr, center, amplitude, hwhm):
    """Return a zero-background Lorentzian peak."""

    return amplitude / (1.0 + ((xarr - center) / hwhm) ** 2)


def _velocity_dho(xarr, center, amplitude, hwhm):
    """Return a peak-normalized velocity-spectrum DHO line shape.

    ``hwhm`` is half of the damping linewidth in ordinary-frequency units, so
    its weak-damping FWHM is ``2*hwhm``.  The current pySED conversion is
    therefore ``tau = 1 / (2*pi*hwhm)`` for both line shapes.
    """

    xarr = np.asarray(xarr, dtype=float)
    damping_width = 2.0 * hwhm
    numerator = amplitude * damping_width**2 * xarr**2
    denominator = (
        (xarr**2 - center**2) ** 2 + damping_width**2 * xarr**2
    )
    result = np.zeros_like(xarr)
    np.divide(
        numerator,
        denominator,
        out=result,
        where=denominator > 0.0,
    )
    result[np.isclose(xarr, center) & np.isclose(denominator, 0.0)] = (
        amplitude
    )
    return result


def _line_shape_function(name):
    """Return the requested single-peak line shape."""

    return {"lorentz": _lorentzian, "dho": _velocity_dho}[name]


def _half_maximum_frequencies(fitting_function, center, hwhm):
    """Return the two half-maximum frequencies of one fitted peak."""

    if fitting_function == "lorentz":
        return center - hwhm, center + hwhm
    if fitting_function == "dho":
        damping_width = 2.0 * hwhm
        root = np.sqrt(4.0 * center**2 + damping_width**2)
        return (
            0.5 * (root - damping_width),
            0.5 * (root + damping_width),
        )
    raise ValueError(f"unsupported fitting function '{fitting_function}'")


def _aicc(residuals, num_parameters):
    """Return the small-sample Akaike information criterion."""

    residuals = np.asarray(residuals, dtype=float)
    num_samples = residuals.size
    if num_samples <= num_parameters + 1:
        return np.inf

    rss = float(np.dot(residuals, residuals))
    mean_square = max(rss / num_samples, np.finfo(float).tiny)
    aic = num_samples * np.log(mean_square) + 2.0 * num_parameters
    correction = (
        2.0
        * num_parameters
        * (num_parameters + 1)
        / (num_samples - num_parameters - 1)
    )
    return aic + correction


def _peak_fit_interval(left_base, right_base, peak, size, modulate_factor):
    """Return the local valley-bounded interval for one peak.

    The interval follows the original pySED slicing rule:
    ``[left_base + modulate_factor, right_base - modulate_factor)``.
    The right endpoint is exclusive, as in the former NumPy slice. A minimal
    local expansion is used only if that range has too few samples for a
    three-parameter fit.
    """

    start = max(0, int(left_base) + modulate_factor)
    end = min(size - 1, int(right_base) - modulate_factor - 1)
    if start >= peak or end <= peak:
        start = max(0, int(peak) - 2)
        end = min(size - 1, int(peak) + 2)

    while end - start + 1 < 5 and (start > 0 or end < size - 1):
        if start > 0:
            start -= 1
        if end - start + 1 >= 5:
            break
        if end < size - 1:
            end += 1
    return start, end


def _local_peak_bases(
    valley_signal,
    peaks,
    peak_number,
):
    """Return the complete smoothed valley-to-valley basin of one peak."""

    valley_signal = np.asarray(valley_signal, dtype=float)
    peaks = np.asarray(peaks, dtype=int)
    peak = int(peaks[peak_number])

    if peak_number > 0:
        previous_peak = int(peaks[peak_number - 1])
        left_base = previous_peak + int(
            np.argmin(valley_signal[previous_peak : peak + 1])
        )
    else:
        left_base = int(np.argmin(valley_signal[: peak + 1]))

    if peak_number + 1 < peaks.size:
        next_peak = int(peaks[peak_number + 1])
        right_base = peak + int(
            np.argmin(valley_signal[peak : next_peak + 1])
        )
    else:
        right_base = peak + int(np.argmin(valley_signal[peak:]))

    return left_base, right_base


def _fit_single_peak_candidate(
    frequency,
    sed,
    peak_index,
    start,
    end,
    initial_hwhm,
    peak_max_hwhm,
    fitting_function,
):
    """Fit one zero-background peak with one prescribed line shape."""

    if fitting_function not in _FITTING_FUNCTIONS:
        raise ValueError(f"unsupported fitting function '{fitting_function}'")

    xarr = np.asarray(frequency[start : end + 1], dtype=float)
    values = np.asarray(sed[start : end + 1], dtype=float)
    if xarr.size < 5:
        raise ValueError("the fitting interval has fewer than five samples")
    spacing = np.diff(xarr)
    positive_spacing = spacing[spacing > 0]
    if not positive_spacing.size:
        raise ValueError("the fitting frequency axis must be increasing")
    frequency_step = float(np.median(positive_spacing))
    center = float(frequency[peak_index])
    if fitting_function == "dho" and center <= frequency_step:
        raise ValueError(
            "velocity DHO is not identifiable for a zero-frequency peak"
        )

    minimum_hwhm = max(frequency_step * 1.0e-6, np.finfo(float).eps)
    maximum_hwhm = max(float(peak_max_hwhm), 10.0 * minimum_hwhm)
    center_lower = float(frequency[max(0, peak_index - 1)])
    center_upper = float(frequency[min(len(frequency) - 1, peak_index + 1)])
    if fitting_function == "dho":
        center_lower = max(0.0, center_lower)
    if center_upper <= center_lower:
        raise ValueError("the peak center cannot be bounded")

    observed_amplitude = max(
        float(sed[peak_index]),
        np.finfo(float).tiny,
    )
    # Fit a dimensionless spectrum so SciPy's convergence tests do not depend
    # on the absolute SED unit scale.  Real SED values are commonly around
    # 1e-6 eV/THz; fitting those values directly can make ``curve_fit`` accept
    # the initial center and HWHM without taking an optimization step.
    scaled_values = values / observed_amplitude
    hwhm_guess = float(
        np.clip(initial_hwhm, 10.0 * minimum_hwhm, 0.95 * maximum_hwhm)
    )
    model = _line_shape_function(fitting_function)
    scaled_optimal, scaled_covariance = curve_fit(
        model,
        xarr,
        scaled_values,
        p0=np.array([center, 1.0, hwhm_guess]),
        bounds=(
            np.array([center_lower, 1.0, minimum_hwhm]),
            np.array(
                [center_upper, 2.0, maximum_hwhm]
            ),
        ),
        maxfev=100000,
    )
    parameter_scale = np.array([1.0, observed_amplitude, 1.0])
    optimal = scaled_optimal * parameter_scale
    covariance = scaled_covariance * np.outer(
        parameter_scale,
        parameter_scale,
    )
    predicted = model(xarr, *optimal)
    residuals = values - predicted
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    hwhm = optimal[2]
    half_maximum = _half_maximum_frequencies(
        fitting_function,
        optimal[0],
        hwhm,
    )
    half_bin = 0.5 * frequency_step
    incomplete_peak_shape = (
        half_maximum[0] < xarr[0] - half_bin
        or half_maximum[1] > xarr[-1] + half_bin
    )
    return {
        "fitting_function": fitting_function,
        "frequency": xarr,
        "observed": values,
        "predicted": predicted,
        "component_curves": np.asarray([predicted]),
        "peak_parameters": optimal[np.newaxis, :],
        "peak_errors": errors[np.newaxis, :],
        "upper_bound_hits": np.asarray([hwhm >= 0.98 * maximum_hwhm]),
        "unresolved_widths": np.asarray([hwhm < 0.5 * frequency_step]),
        "incomplete_peak_shapes": np.asarray([incomplete_peak_shape]),
        "half_maximum_frequencies": np.asarray(half_maximum),
        "aicc": _aicc(residuals, 3),
        "rss": float(np.dot(residuals, residuals)),
        "fit_start": float(xarr[0]),
        "fit_end": float(xarr[-1]),
        "num_points": xarr.size,
    }


def _fit_single_peak(
    frequency,
    sed,
    peak_index,
    start,
    end,
    initial_hwhm=0.001,
    peak_max_hwhm=1.0e6,
    fitting_function="auto",
):
    """Fit one peak, selecting Lorentz/DHO by AICc when requested."""

    if fitting_function not in {*_FITTING_FUNCTIONS, "auto"}:
        raise ValueError(f"unsupported fitting function '{fitting_function}'")
    requested_functions = (
        _FITTING_FUNCTIONS if fitting_function == "auto" else (fitting_function,)
    )
    candidates = []
    errors = []
    for candidate_function in requested_functions:
        try:
            candidates.append(
                _fit_single_peak_candidate(
                    frequency,
                    sed,
                    peak_index,
                    start,
                    end,
                    initial_hwhm,
                    peak_max_hwhm,
                    candidate_function,
                )
            )
        except (RuntimeError, TypeError, ValueError) as error:
            errors.append(f"{candidate_function}: {error}")

    if not candidates:
        raise RuntimeError("; ".join(errors) or "all peak fits failed")
    if fitting_function != "auto":
        return candidates[0]

    lorentz_candidates = [
        candidate
        for candidate in candidates
        if candidate["fitting_function"] == "lorentz"
    ]
    dho_candidates = [
        candidate
        for candidate in candidates
        if candidate["fitting_function"] == "dho"
    ]
    # A DHO becomes numerically indistinguishable from a Lorentzian in the
    # weak-damping limit. Prefer the simpler conventional interpretation.
    if lorentz_candidates and dho_candidates:
        dho_ratio = (
            dho_candidates[0]["peak_parameters"][0, 2]
            / dho_candidates[0]["peak_parameters"][0, 0]
        )
        if dho_ratio < _DHO_MIN_DISTINGUISHABLE_DAMPING_RATIO:
            return lorentz_candidates[0]

    best_aicc = min(candidate["aicc"] for candidate in candidates)
    statistically_equivalent = [
        candidate for candidate in candidates
        if candidate["aicc"] <= best_aicc + _AUTO_AICC_EQUIVALENCE
    ]
    return min(
        statistically_equivalent,
        key=lambda candidate: 0
        if candidate["fitting_function"] == "lorentz"
        else 1,
    )


class lorentz:
    """Detect and independently fit zero-background SED peaks."""

    def __init__(self, data, params):
        self.q_index = params.qpoint_slice_index
        self.lorentz_fit_freq_min = getattr(params, "lorentz_fit_freq_min", None)
        self.lorentz_fit_freq_max = getattr(params, "lorentz_fit_freq_max", None)
        self.fitting_function = getattr(params, "fitting_function", "auto")

        selected_frequency, selected_sed = _select_lorentz_frequency_range(
            data.freq_fft,
            data.sed_avg[:, self.q_index],
            freq_min=self.lorentz_fit_freq_min,
            freq_max=self.lorentz_fit_freq_max,
        )
        if (
            self.lorentz_fit_freq_min is not None
            or self.lorentz_fit_freq_max is not None
        ):
            print(
                "      Frequency range      : "
                f"{selected_frequency.min():.3f}-"
                f"{selected_frequency.max():.3f} THz\n"
            )

        # Mirror only a spectrum that starts at zero.  This keeps a possible
        # acoustic endpoint maximum detectable while never fitting it twice.
        mirror_points = min(5, selected_sed.size - 1)
        starts_at_zero = np.isclose(selected_frequency[0], 0.0)
        if mirror_points and starts_at_zero:
            self.sed = np.concatenate(
                [selected_sed[1 : mirror_points + 1][::-1], selected_sed]
            )
            self.thz = np.concatenate(
                [
                    -selected_frequency[1 : mirror_points + 1][::-1],
                    selected_frequency,
                ]
            )
        else:
            self.sed = np.asarray(selected_sed, dtype=float)
            self.thz = np.asarray(selected_frequency, dtype=float)

        peaks, peak_significance = _find_sed_peaks(
            self.sed,
            peak_min_significance=params.peak_min_significance,
        )
        keep = self.thz[peaks] >= 0.0
        peaks = peaks[keep]
        peak_significance = peak_significance[keep]
        valley_signal = _smoothed_log_sed(self.sed)

        print("  >  Peak detection")
        print(f"      Peaks found          : {len(peaks)}")
        print(
            "      Minimum significance : "
            f"{params.peak_min_significance:g}"
        )
        if self.fitting_function == "auto":
            print("      Line shape           : automatic Lorentz or DHO")
        else:
            print(f"      Line shape           : {self.fitting_function}")
        fit_figure = Plot_SED.resolve_slice_output_path(params, lorentz=True)
        print(f"      Fit figure           : {fit_figure}")

        self.fit_clusters = []
        for peak_number, peak in enumerate(peaks, start=1):
            left_base, right_base = _local_peak_bases(
                valley_signal,
                peaks,
                peak_number - 1,
            )
            start, end = _peak_fit_interval(
                left_base,
                right_base,
                peak,
                self.sed.size,
                params.modulate_factor,
            )
            try:
                result = _fit_single_peak(
                    self.thz,
                    self.sed,
                    peak,
                    start,
                    end,
                    initial_hwhm=params.initial_guess_hwhm,
                    peak_max_hwhm=params.peak_max_hwhm,
                    fitting_function=self.fitting_function,
                )
            except RuntimeError as error:
                print(
                    "  WARNING: spectrum fit failed for peak at "
                    f"{self.thz[peak]:.6f} THz: {error}"
                )
                continue

            result["peak_significance"] = float(
                peak_significance[peak_number - 1]
            )
            result["peak_number"] = peak_number
            self.fit_clusters.append(result)

        print("\n  >  Independent peak fits")
        range_width = 19
        print(
            f"      {'No.':>3}   {'Fit range (THz)':^{range_width}}   "
            f"{'Data points':>11}"
        )
        print(
            f"      {'---':>3}   {'-' * range_width}   "
            f"{'-' * 11}"
        )
        for result in self.fit_clusters:
            fit_range = (
                f"{result['fit_start']:8.4f} - "
                f"{result['fit_end']:<8.4f}"
            )
            print(
                f"      {result['peak_number']:3d}   "
                f"{fit_range}   "
                f"{result['num_points']:11d}"
            )
        upper_bound_peaks = [
            str(result["peak_number"])
            for result in self.fit_clusters
            if result["upper_bound_hits"][0]
        ]
        incomplete_peaks = [
            str(result["peak_number"])
            for result in self.fit_clusters
            if result["incomplete_peak_shapes"][0]
        ]
        if upper_bound_peaks:
            print(
                "  Note: fitted HWHM reached peak_max_hwhm for peak(s): "
                + ", ".join(upper_bound_peaks)
            )
        if incomplete_peaks:
            peak_list = ", ".join(incomplete_peaks)
            if len(incomplete_peaks) == 1:
                print(
                    f"  Note: peak {peak_list} is incomplete; HWHM is "
                    "model-extrapolated."
                )
                print(
                    "        Lifetime is qualitative. Adjust "
                    "peak_min_significance if needed."
                )
            else:
                print(
                    f"  Note: peaks {peak_list} are incomplete; HWHMs are "
                    "model-extrapolated."
                )
                print(
                    "        Lifetimes are qualitative. Adjust "
                    "peak_min_significance if needed."
                )

        if self.fit_clusters:
            self.popt = np.vstack(
                [result["peak_parameters"] for result in self.fit_clusters]
            )
            self.pcov = np.vstack(
                [result["peak_errors"] for result in self.fit_clusters]
            )
            self.fit_models = np.asarray(
                [result["fitting_function"] for result in self.fit_clusters],
                dtype=object,
            )
            self.upper_bound_hits = np.asarray(
                [result["upper_bound_hits"][0] for result in self.fit_clusters],
                dtype=bool,
            )
            self.unresolved_widths = np.asarray(
                [result["unresolved_widths"][0] for result in self.fit_clusters],
                dtype=bool,
            )
            self.incomplete_peak_shapes = np.asarray(
                [
                    result["incomplete_peak_shapes"][0]
                    for result in self.fit_clusters
                ],
                dtype=bool,
            )
            self.peak_significance = np.asarray(
                [result["peak_significance"] for result in self.fit_clusters],
                dtype=float,
            )
            order = np.argsort(self.popt[:, 0])
            self.popt = self.popt[order]
            self.pcov = self.pcov[order]
            self.fit_models = self.fit_models[order]
            self.upper_bound_hits = self.upper_bound_hits[order]
            self.unresolved_widths = self.unresolved_widths[order]
            self.incomplete_peak_shapes = self.incomplete_peak_shapes[order]
            self.peak_significance = self.peak_significance[order]
            self.fit_clusters = [self.fit_clusters[index] for index in order]
        else:
            self.popt = np.empty((0, 3), dtype=float)
            self.pcov = np.empty((0, 3), dtype=float)
            self.fit_models = np.empty(0, dtype=object)
            self.upper_bound_hits = np.empty(0, dtype=bool)
            self.unresolved_widths = np.empty(0, dtype=bool)
            self.incomplete_peak_shapes = np.empty(0, dtype=bool)
            self.peak_significance = np.empty(0, dtype=float)

        params.popt = self.popt
        params.pcov = self.pcov
        params.fit_models = self.fit_models
        params.fit_peak_significance = self.peak_significance
        params.lorentz_fit_clusters = self.fit_clusters
        params.plot_lorentz = True

        lifetime_file = FileIO.write_phonon_lifetime(self, params)

        lifetimes = FileIO.hwhm_to_lifetime_ps(self.popt[:, 2])
        print("\n  >  Spectral fitting results")
        print("      Lifetime definition  : tau_SED = 1 / (2*pi*HWHM)")
        print(
            "      Frequency (THz)      HWHM (THz)      Tau_SED (ps)"
            "        Model"
        )
        print(
            "      ---------------      ----------      -----------"
            "        -------"
        )
        for frequency, hwhm, lifetime, model in zip(
            self.popt[:, 0],
            self.popt[:, 2],
            lifetimes,
            self.fit_models,
        ):
            if np.isfinite(lifetime):
                print(
                    f"      {frequency:15.6f}      {hwhm:10.6f}      "
                    f"{lifetime:11.8f}        {model}"
                )
        overdamped_dho = (
            (self.fit_models == "dho")
            & (self.popt[:, 2] >= self.popt[:, 0])
        )
        if np.any(overdamped_dho):
            print(
                "  Note: overdamped DHO lifetimes are included for "
                "qualitative comparison only."
            )
        print(f"\n  [OK] Lifetime data written: {lifetime_file}")
        Plot_SED.plot_slice(data, params)

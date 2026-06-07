"""Quantum detailed-balance corrections for classical scattering spectra."""

import numpy as np


PLANCK_MEV_THZ = 4.135667696
"""Planck constant expressed as meV per THz."""

BOLTZMANN_MEV_PER_K = 0.08617333262145
"""Boltzmann constant expressed as meV/K."""

WAVENUMBER_MEV = 0.1239841984332003
"""Energy of one inverse centimeter expressed in meV."""


def energy_axis_to_mev(axis, unit="thz"):
    """Convert an energy/frequency axis to signed meV values."""

    values = np.asarray(axis, dtype=float)
    unit_key = str(unit).lower().replace("_", "-")
    if unit_key in ("mev", "millielectronvolt", "millielectronvolts"):
        return values
    if unit_key in ("ev", "electronvolt", "electronvolts"):
        return values * 1000.0
    if unit_key in ("thz", "frequency-thz", "hz-thz"):
        return values * PLANCK_MEV_THZ
    if unit_key in ("cm-1", "cm^-1", "wavenumber", "wavenumbers"):
        return values * WAVENUMBER_MEV
    raise ValueError("unit must be 'thz', 'meV', 'eV', or 'cm-1'")


def energy_axis_from_mev(energy_mev, unit="thz"):
    """Convert signed meV values to the requested unit."""

    values = np.asarray(energy_mev, dtype=float)
    unit_key = str(unit).lower().replace("_", "-")
    if unit_key in ("mev", "millielectronvolt", "millielectronvolts"):
        return values
    if unit_key in ("ev", "electronvolt", "electronvolts"):
        return values / 1000.0
    if unit_key in ("thz", "frequency-thz", "hz-thz"):
        return values / PLANCK_MEV_THZ
    if unit_key in ("cm-1", "cm^-1", "wavenumber", "wavenumbers"):
        return values / WAVENUMBER_MEV
    raise ValueError("unit must be 'thz', 'meV', 'eV', or 'cm-1'")


def convert_energy_axis(axis, source_unit="thz", target_unit="meV"):
    """Convert an axis between supported frequency/energy units."""

    return energy_axis_from_mev(
        energy_axis_to_mev(axis, unit=source_unit),
        unit=target_unit,
    )


def classical_to_quantum_factor(axis, temperature, unit="thz"):
    """Return the harmonic classical-to-quantum correction factor.

    The factor is

    ``x / (1 - exp(-x))``, where ``x = E / (k_B T)``.

    For negative energy transfers the same signed expression obeys detailed
    balance.  The zero-energy limit is evaluated with the Taylor expansion.
    """

    if temperature is None:
        return np.ones_like(np.asarray(axis, dtype=float))
    if temperature <= 0:
        raise ValueError("temperature must be positive")

    energy_mev = energy_axis_to_mev(axis, unit=unit)
    x = energy_mev / (BOLTZMANN_MEV_PER_K * float(temperature))
    factor = np.empty_like(x, dtype=float)
    small = np.abs(x) < 1e-8
    factor[small] = 1.0 + 0.5 * x[small] + (x[small] ** 2) / 12.0
    with np.errstate(over="ignore", invalid="ignore"):
        factor[~small] = x[~small] / (-np.expm1(-x[~small]))
    return factor


def bose_population(axis, temperature, unit="thz"):
    """Return the Bose-Einstein population for positive energy magnitudes."""

    if temperature <= 0:
        raise ValueError("temperature must be positive")
    energy_mev = np.abs(energy_axis_to_mev(axis, unit=unit))
    x = energy_mev / (BOLTZMANN_MEV_PER_K * float(temperature))
    out = np.empty_like(x, dtype=float)
    zero = x == 0
    out[zero] = np.inf
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        out[~zero] = 1.0 / np.expm1(x[~zero])
    return out


def apply_quantum_correction(intensity, axis, temperature, unit="thz", axis_index=-1):
    """Multiply a spectrum or map by the harmonic quantum correction factor."""

    arr = np.asarray(intensity, dtype=float)
    factor = classical_to_quantum_factor(axis, temperature, unit=unit)
    if factor.ndim != 1:
        raise ValueError("axis must be one-dimensional")
    axis_index = int(axis_index)
    if axis_index < 0:
        axis_index += arr.ndim
    if axis_index < 0 or axis_index >= arr.ndim:
        raise ValueError("axis_index is outside intensity dimensions")
    if arr.shape[axis_index] != factor.size:
        raise ValueError("axis length must match intensity dimension")
    shape = [1] * arr.ndim
    shape[axis_index] = factor.size
    return arr * factor.reshape(shape)

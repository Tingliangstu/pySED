"""Dynamic structure factor calculations from MD trajectories.

This module implements the density-correlation route used for neutron and
X-ray inelastic scattering.  The computation is organized in terms of

    rho(Q, t) = sum_i w_i(Q) exp(i Q dot r_i(t))

and estimates the coherent spectrum from the time Fourier transform of rho.
This is equivalent to evaluating the intermediate scattering function and then
Fourier transforming it, but avoids storing F(Q, t) explicitly.
"""

from dataclasses import dataclass

import numpy as np

from pySED.probe_weights import ATOMIC_NUMBERS, atomic_number_weights
from pySED.probe_weights import xray_form_factor_weights as _xray_form_factor_weights
from pySED.scattering_kernel import (
    as_blocks,
    circular_autocorrelation,
    density_amplitude,
    density_amplitude_from_phase,
    fft_power,
    frequency_grid,
    phase_factors,
    resolve_backend,
    standard_error,
    to_numpy,
)


NEUTRON_COHERENT_LENGTHS_FM = {
    "H": -3.7390,
    "D": 6.671,
    "B": 5.30,
    "C": 6.6460,
    "N": 9.36,
    "O": 5.803,
    "Na": 3.63,
    "Al": 3.449,
    "Si": 4.1491,
    "P": 5.13,
    "S": 2.847,
    "Mo": 6.715,
    "Cu": 7.718,
    "Se": 7.970,
    "Zr": 7.16,
    "W": 4.86,
}

NEUTRON_INCOHERENT_XS_BARN = {
    "H": 80.26,
    "D": 2.05,
    "B": 1.7,
    "C": 0.001,
    "N": 0.5,
    "O": 0.0008,
    "Na": 1.62,
    "Al": 0.0082,
    "Si": 0.004,
    "P": 0.005,
    "S": 0.007,
    "Mo": 0.04,
    "Cu": 0.55,
    "Se": 0.32,
    "Zr": 0.02,
    "W": 1.63,
}

@dataclass
class DSFResult:
    qpoints_cartesian: np.ndarray
    frequencies_thz: np.ndarray
    coherent: np.ndarray
    incoherent: np.ndarray
    total: np.ndarray
    elastic_coherent: np.ndarray
    elastic_incoherent: np.ndarray
    num_blocks: int
    metadata: dict
    coherent_sem: np.ndarray = None
    incoherent_sem: np.ndarray = None
    total_sem: np.ndarray = None
    elastic_coherent_sem: np.ndarray = None
    elastic_incoherent_sem: np.ndarray = None


@dataclass
class CurrentCorrelationResult:
    qpoints_cartesian: np.ndarray
    frequencies_thz: np.ndarray
    longitudinal: np.ndarray
    transverse: np.ndarray
    num_blocks: int
    longitudinal_sem: np.ndarray = None
    transverse_sem: np.ndarray = None
    metadata: dict = None


@dataclass
class PartialDSFResult:
    qpoints_cartesian: np.ndarray
    frequencies_thz: np.ndarray
    species: tuple
    partial: np.ndarray
    partial_sem: np.ndarray
    weighted_total: np.ndarray
    weighted_total_sem: np.ndarray
    num_blocks: int
    metadata: dict


def _weights_from_atom_types(atom_types, table, missing="raise"):
    if atom_types is None:
        return None
    values = []
    for atom in atom_types:
        key = str(atom)
        if key not in table:
            if missing == "zero":
                values.append(0.0)
                continue
            raise KeyError("missing scattering data for atom type %s" % key)
        values.append(float(table[key]))
    return np.asarray(values, dtype=float)


def neutron_weights(atom_types, include_incoherent=True):
    coherent = _weights_from_atom_types(atom_types, NEUTRON_COHERENT_LENGTHS_FM)
    incoherent = None
    if include_incoherent:
        incoherent = _weights_from_atom_types(atom_types, NEUTRON_INCOHERENT_XS_BARN)
    return coherent, incoherent


def xray_atomic_number_weights(atom_types, qpoints_cartesian=None):
    """Return an approximate X-ray weight using atomic number.

    This is a safe default when no tabulated form-factor data are supplied.  For
    publication-quality X-ray comparison, pass q-dependent form factors through
    ``coherent_weights``.
    """

    return atomic_number_weights(atom_types, qpoints_cartesian=qpoints_cartesian)


def xray_form_factor_weights(
    atom_types,
    qpoints_cartesian=None,
    coefficients_table=None,
    missing="atomic_number",
):
    """Return q-dependent X-ray form-factor weights.

    This thin wrapper exposes :mod:`pySED.probe_weights` through the historical
    DSF namespace.
    """

    return _xray_form_factor_weights(
        atom_types,
        qpoints_cartesian=qpoints_cartesian,
        coefficients_table=coefficients_table,
        missing=missing,
    )


def resolve_scattering_weights(
    atom_types=None,
    qpoints_cartesian=None,
    experiment="neutron",
    coherent_weights=None,
    incoherent_cross_sections=None,
    xray_form_factor_table=None,
    xray_missing="atomic_number",
):
    qpoints = np.atleast_2d(qpoints_cartesian) if qpoints_cartesian is not None else None

    if coherent_weights is not None:
        coherent = np.asarray(coherent_weights, dtype=float)
    elif experiment == "neutron":
        coherent, _ = neutron_weights(atom_types, include_incoherent=False)
    elif experiment == "xray":
        coherent = xray_form_factor_weights(
            atom_types,
            qpoints_cartesian=qpoints,
            coefficients_table=xray_form_factor_table,
            missing=xray_missing,
        )
    elif experiment == "custom":
        coherent = None
    else:
        raise ValueError("experiment must be 'neutron', 'xray', or 'custom'")

    if incoherent_cross_sections is not None:
        incoherent = np.asarray(incoherent_cross_sections, dtype=float)
    elif experiment == "neutron":
        _, incoherent = neutron_weights(atom_types, include_incoherent=True)
    else:
        incoherent = None

    return coherent, incoherent


def _weights_for_q(weights, q_index, num_atoms):
    if weights is None:
        return np.ones(num_atoms, dtype=float)
    weights = np.asarray(weights, dtype=float)
    if weights.ndim == 1:
        if weights.shape[0] != num_atoms:
            raise ValueError("weights length must match number of atoms")
        return weights
    if weights.ndim == 2:
        if weights.shape[1] != num_atoms:
            raise ValueError("q-dependent weights must have shape (num_qpoints, num_atoms)")
        return weights[q_index]
    raise ValueError("weights must be one- or two-dimensional")


def _species_groups(atom_types):
    if atom_types is None:
        raise ValueError("atom_types must be supplied for partial DSF")
    labels = tuple(dict.fromkeys(str(atom) for atom in atom_types))
    atom_types = np.asarray([str(atom) for atom in atom_types])
    groups = [np.flatnonzero(atom_types == label) for label in labels]
    return labels, groups


def _species_weights_for_q(species_weights, species, q_index):
    n_species = len(species)
    if species_weights is None:
        return np.ones(n_species, dtype=float)
    if isinstance(species_weights, dict):
        return np.asarray([species_weights[label] for label in species])
    weights = np.asarray(species_weights)
    if weights.ndim == 1:
        if weights.shape[0] != n_species:
            raise ValueError("species_weights length must match number of species")
        return weights
    if weights.ndim == 2:
        if weights.shape[1] != n_species:
            raise ValueError("q-dependent species_weights must have shape (num_qpoints, num_species)")
        return weights[q_index]
    raise ValueError("species_weights must be a dict, 1D array, or 2D array")


def compute_dsf(
    positions,
    qpoints_cartesian,
    dt,
    atom_types=None,
    experiment="neutron",
    coherent_weights=None,
    incoherent_cross_sections=None,
    xray_form_factor_table=None,
    xray_missing="atomic_number",
    num_blocks=1,
    positive_only=True,
    backend="cpu",
):
    """Compute coherent/incoherent dynamic structure factors.

    Parameters
    ----------
    positions
        Array with shape ``(num_frames, num_atoms, 3)`` in Angstrom.
    qpoints_cartesian
        Q points in inverse Angstrom, including the ``2*pi`` convention.
    dt
        Sampling interval in seconds.
    backend
        ``"cpu"``/``"numpy"`` by default.  Use ``"cupy"`` or ``"gpu"`` only
        when CuPy is installed for the local CUDA runtime.
    """

    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 3 or positions.shape[2] != 3:
        raise ValueError("positions must have shape (num_frames, num_atoms, 3)")
    qpoints = np.atleast_2d(np.asarray(qpoints_cartesian, dtype=float))
    if qpoints.shape[1] != 3:
        raise ValueError("qpoints_cartesian must have shape (num_qpoints, 3)")
    kernel = resolve_backend(backend)
    xp = kernel.xp

    n_frames, n_atoms, _ = positions.shape
    blocks = as_blocks(n_frames, num_blocks)
    block_size = blocks[0][1] - blocks[0][0]
    grid = frequency_grid(block_size, dt, positive_only=positive_only, backend=kernel)
    freq_slice = grid.freq_slice
    freq_out = grid.frequencies_thz
    user_coherent_weights = coherent_weights is not None

    coherent_weights, incoherent_cross_sections = resolve_scattering_weights(
        atom_types=atom_types,
        qpoints_cartesian=qpoints,
        experiment=experiment,
        coherent_weights=coherent_weights,
        incoherent_cross_sections=incoherent_cross_sections,
        xray_form_factor_table=xray_form_factor_table,
        xray_missing=xray_missing,
    )

    n_q = qpoints.shape[0]
    n_f = len(freq_out)
    coherent_blocks = np.zeros((num_blocks, n_q, n_f), dtype=float)
    incoherent_blocks = np.zeros((num_blocks, n_q, n_f), dtype=float)
    elastic_coherent_blocks = np.zeros((num_blocks, n_q), dtype=float)
    elastic_incoherent_blocks = np.zeros((num_blocks, n_q), dtype=float)

    calc_coherent = coherent_weights is not None or experiment in ("neutron", "xray", "custom")
    calc_incoherent = incoherent_cross_sections is not None

    block_norm = float(block_size * n_atoms)
    for block_index, (start, stop) in enumerate(blocks):
        block = positions[start:stop]
        for q_index, q_vec in enumerate(qpoints):
            phase = phase_factors(block, q_vec, backend=kernel)

            if calc_coherent:
                w = _weights_for_q(coherent_weights, q_index, n_atoms)
                rho = density_amplitude_from_phase(phase, w, backend=kernel)
                elastic_coherent_blocks[block_index, q_index] = (
                    to_numpy(xp.abs(xp.mean(rho)) ** 2, kernel) / block_norm
                )
                spectrum = fft_power(rho, norm="forward", backend=kernel)
                coherent_blocks[block_index, q_index] = (
                    to_numpy(spectrum[freq_slice].real, kernel) / block_norm
                )

            if calc_incoherent:
                sigma = _weights_for_q(incoherent_cross_sections, q_index, n_atoms)
                sigma_backend = xp.asarray(sigma)
                atom_fft = kernel.fft(phase.T, axis=1, norm="forward")
                elastic_incoherent_blocks[block_index, q_index] = (
                    to_numpy(
                        xp.sum(sigma_backend * xp.abs(xp.mean(phase.T, axis=1)) ** 2),
                        kernel,
                    )
                    / block_norm
                )
                spectrum = xp.sum(
                    sigma_backend.reshape(-1, 1) * xp.abs(atom_fft) ** 2,
                    axis=0,
                )
                incoherent_blocks[block_index, q_index] = (
                    to_numpy(spectrum[freq_slice].real, kernel) / block_norm
                )

    coherent = np.mean(coherent_blocks, axis=0)
    incoherent = np.mean(incoherent_blocks, axis=0)
    elastic_coherent = np.mean(elastic_coherent_blocks, axis=0)
    elastic_incoherent = np.mean(elastic_incoherent_blocks, axis=0)
    coherent_sem = standard_error(coherent_blocks)
    incoherent_sem = standard_error(incoherent_blocks)
    elastic_coherent_sem = standard_error(elastic_coherent_blocks)
    elastic_incoherent_sem = standard_error(elastic_incoherent_blocks)
    total_blocks = coherent_blocks + incoherent_blocks
    total = coherent + incoherent
    total_sem = standard_error(total_blocks)

    return DSFResult(
        qpoints_cartesian=qpoints,
        frequencies_thz=freq_out,
        coherent=coherent,
        incoherent=incoherent,
        total=total,
        elastic_coherent=elastic_coherent,
        elastic_incoherent=elastic_incoherent,
        num_blocks=num_blocks,
        metadata={
            "experiment": experiment,
            "dt_seconds": float(dt),
            "positive_only": bool(positive_only),
            "backend": kernel.name,
            "xray_weights": (
                "user supplied coherent_weights"
                if experiment == "xray" and user_coherent_weights
                else "Cromer-Mann form factors with atomic-number fallback"
                if experiment == "xray"
                else None
            ),
            "normalization": "block spectra divided by block_size * num_atoms, then averaged",
            "uncertainty": "standard error of the block-averaged spectra",
        },
        coherent_sem=coherent_sem,
        incoherent_sem=incoherent_sem,
        total_sem=total_sem,
        elastic_coherent_sem=elastic_coherent_sem,
        elastic_incoherent_sem=elastic_incoherent_sem,
    )


def compute_partial_dsf(
    positions,
    qpoints_cartesian,
    dt,
    atom_types,
    species_weights=None,
    num_blocks=1,
    positive_only=True,
    backend="cpu",
):
    """Compute species-pair coherent partial DSF matrices.

    The returned ``partial`` array has shape
    ``(num_qpoints, num_species, num_species, num_frequencies)`` and contains
    the block-averaged cross spectra

    ``FFT[rho_A(Q,t)] * conj(FFT[rho_B(Q,t)))``.

    ``weighted_total`` is the coherent spectrum reconstructed from the partials
    using the supplied species weights.
    """

    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 3 or positions.shape[2] != 3:
        raise ValueError("positions must have shape (num_frames, num_atoms, 3)")
    qpoints = np.atleast_2d(np.asarray(qpoints_cartesian, dtype=float))
    if qpoints.shape[1] != 3:
        raise ValueError("qpoints_cartesian must have shape (num_qpoints, 3)")

    n_frames, n_atoms, _ = positions.shape
    atom_types = np.asarray(atom_types)
    if atom_types.shape[0] != n_atoms:
        raise ValueError("atom_types length must match number of atoms")

    kernel = resolve_backend(backend)
    xp = kernel.xp
    species, groups = _species_groups(atom_types)
    n_species = len(species)
    blocks = as_blocks(n_frames, num_blocks)
    block_size = blocks[0][1] - blocks[0][0]
    grid = frequency_grid(block_size, dt, positive_only=positive_only, backend=kernel)
    freq_slice = grid.freq_slice
    freq_out = grid.frequencies_thz

    n_q = qpoints.shape[0]
    n_f = len(freq_out)
    partial_blocks = np.zeros((num_blocks, n_q, n_species, n_species, n_f), dtype=complex)
    weighted_total_blocks = np.zeros((num_blocks, n_q, n_f), dtype=float)
    block_norm = float(block_size * n_atoms)

    for block_index, (start, stop) in enumerate(blocks):
        block = positions[start:stop]
        for q_index, q_vec in enumerate(qpoints):
            phase = phase_factors(block, q_vec, backend=kernel)
            rho_species = xp.zeros((block_size, n_species), dtype=complex)
            for species_index, atom_ids in enumerate(groups):
                rho_species[:, species_index] = xp.sum(
                    phase[:, xp.asarray(atom_ids)],
                    axis=1,
                )

            rho_fft = kernel.fft(rho_species, axis=0, norm="forward")[freq_slice]
            partial = xp.einsum("fa,fb->abf", rho_fft, xp.conjugate(rho_fft)) / block_norm
            partial_np = to_numpy(partial, kernel)
            partial_blocks[block_index, q_index] = partial_np

            weights = _species_weights_for_q(species_weights, species, q_index)
            weighted = np.einsum("a,abf,b->f", weights, partial_np, np.conjugate(weights))
            weighted_total_blocks[block_index, q_index] = weighted.real

    partial_mean = np.mean(partial_blocks, axis=0)
    partial_sem = np.std(partial_blocks, axis=0, ddof=1) / np.sqrt(num_blocks) if num_blocks > 1 else np.zeros_like(partial_mean, dtype=float)
    weighted_total = np.mean(weighted_total_blocks, axis=0)
    weighted_total_sem = standard_error(weighted_total_blocks)

    return PartialDSFResult(
        qpoints_cartesian=qpoints,
        frequencies_thz=freq_out,
        species=species,
        partial=partial_mean,
        partial_sem=partial_sem,
        weighted_total=weighted_total,
        weighted_total_sem=weighted_total_sem,
        num_blocks=num_blocks,
        metadata={
            "dt_seconds": float(dt),
            "positive_only": bool(positive_only),
            "backend": kernel.name,
            "normalization": "block spectra divided by block_size * num_atoms, then averaged",
            "species_weights": "unit weights" if species_weights is None else "user supplied",
        },
    )


def compute_coherent_dsf_via_correlation(
    positions,
    qpoints_cartesian,
    dt,
    coherent_weights=None,
    num_blocks=1,
    positive_only=True,
    backend="cpu",
):
    """Reference coherent DSF from the density autocorrelation function.

    This routine is intentionally more explicit and slower than
    :func:`compute_dsf`.  It first constructs the finite-block circular
    intermediate scattering function

    ``F(Q, lag) = mean_t rho(Q, t + lag) rho(Q, t).conjugate()``

    and then Fourier transforms ``F`` in time.  It is mainly useful for
    validating that the production estimator based on
    ``abs(FFT_t[rho(Q, t)])**2`` is consistent with the correlation-function
    formulation.
    """

    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 3 or positions.shape[2] != 3:
        raise ValueError("positions must have shape (num_frames, num_atoms, 3)")
    qpoints = np.atleast_2d(np.asarray(qpoints_cartesian, dtype=float))
    if qpoints.shape[1] != 3:
        raise ValueError("qpoints_cartesian must have shape (num_qpoints, 3)")
    kernel = resolve_backend(backend)
    xp = kernel.xp

    n_frames, n_atoms, _ = positions.shape
    blocks = as_blocks(n_frames, num_blocks)
    block_size = blocks[0][1] - blocks[0][0]
    grid = frequency_grid(block_size, dt, positive_only=positive_only, backend=kernel)
    freq_slice = grid.freq_slice
    freq_out = grid.frequencies_thz

    n_q = qpoints.shape[0]
    n_f = len(freq_out)
    coherent_blocks = np.zeros((num_blocks, n_q, n_f), dtype=float)
    elastic_blocks = np.zeros((num_blocks, n_q), dtype=float)
    block_norm = float(block_size * n_atoms)

    for block_index, (start, stop) in enumerate(blocks):
        block = positions[start:stop]
        for q_index, q_vec in enumerate(qpoints):
            w = _weights_for_q(coherent_weights, q_index, n_atoms)
            rho = density_amplitude(block, q_vec, w, backend=kernel)
            elastic_blocks[block_index, q_index] = (
                to_numpy(xp.abs(xp.mean(rho)) ** 2, kernel) / block_norm
            )

            correlation = circular_autocorrelation(rho, backend=kernel)
            spectrum = kernel.fft(correlation, norm="forward").real
            coherent_blocks[block_index, q_index] = (
                to_numpy(spectrum[freq_slice], kernel) / block_norm
            )

    coherent = np.mean(coherent_blocks, axis=0)
    zeros = np.zeros_like(coherent)
    elastic = np.mean(elastic_blocks, axis=0)
    elastic_zeros = np.zeros_like(elastic)
    coherent_sem = standard_error(coherent_blocks)

    return DSFResult(
        qpoints_cartesian=qpoints,
        frequencies_thz=freq_out,
        coherent=coherent,
        incoherent=zeros,
        total=coherent,
        elastic_coherent=elastic,
        elastic_incoherent=elastic_zeros,
        num_blocks=num_blocks,
        metadata={
            "experiment": "custom",
            "dt_seconds": float(dt),
            "positive_only": bool(positive_only),
            "backend": kernel.name,
            "estimator": "circular density autocorrelation followed by time FFT",
            "normalization": "block spectra divided by block_size * num_atoms, then averaged",
            "uncertainty": "standard error of the block-averaged spectra",
        },
        coherent_sem=coherent_sem,
        incoherent_sem=zeros,
        total_sem=coherent_sem,
        elastic_coherent_sem=standard_error(elastic_blocks),
        elastic_incoherent_sem=elastic_zeros,
    )


def compute_current_correlations(
    positions,
    velocities,
    qpoints_cartesian,
    dt,
    weights=None,
    num_blocks=1,
    positive_only=True,
    backend="cpu",
):
    """Compute longitudinal and transverse current spectra.

    Spectra are evaluated independently for each time block, normalized by
    ``block_size * num_atoms``, then averaged.  When ``num_blocks > 1``, the
    result includes standard errors of the block-averaged spectra.
    """

    positions = np.asarray(positions, dtype=float)
    velocities = np.asarray(velocities, dtype=float)
    qpoints = np.atleast_2d(np.asarray(qpoints_cartesian, dtype=float))
    if positions.shape != velocities.shape:
        raise ValueError("positions and velocities must have identical shape")
    if positions.ndim != 3 or positions.shape[2] != 3:
        raise ValueError("positions must have shape (num_frames, num_atoms, 3)")
    kernel = resolve_backend(backend)
    xp = kernel.xp

    n_frames, n_atoms, _ = positions.shape
    blocks = as_blocks(n_frames, num_blocks)
    block_size = blocks[0][1] - blocks[0][0]
    grid = frequency_grid(block_size, dt, positive_only=positive_only, backend=kernel)
    freq_slice = grid.freq_slice
    freq_out = grid.frequencies_thz

    n_q = qpoints.shape[0]
    n_f = len(freq_out)
    longitudinal_blocks = np.zeros((num_blocks, n_q, n_f), dtype=float)
    transverse_blocks = np.zeros_like(longitudinal_blocks)

    block_norm = float(block_size * n_atoms)
    for block_index, (start, stop) in enumerate(blocks):
        pos_block = positions[start:stop]
        vel_block = velocities[start:stop]
        for q_index, q_vec in enumerate(qpoints):
            q_norm = np.linalg.norm(q_vec)
            if q_norm == 0:
                continue
            q_hat = xp.asarray(q_vec / q_norm)
            phase = phase_factors(pos_block, q_vec, backend=kernel)
            w = xp.asarray(_weights_for_q(weights, q_index, n_atoms))
            weighted_phase = phase * w.reshape(1, -1)
            vel_backend = xp.asarray(vel_block, dtype=float)

            v_long = xp.tensordot(vel_backend, q_hat, axes=([2], [0]))
            j_long = xp.sum(weighted_phase * v_long, axis=1)
            long_spec = fft_power(j_long, norm="forward", backend=kernel)
            longitudinal_blocks[block_index, q_index] = (
                to_numpy(long_spec[freq_slice].real, kernel) / block_norm
            )

            v_long_vec = v_long[:, :, None] * q_hat.reshape(1, 1, 3)
            v_trans = vel_backend - v_long_vec
            j_trans = xp.sum(weighted_phase[:, :, None] * v_trans, axis=1)
            trans_spec = xp.sum(fft_power(j_trans, axis=0, norm="forward", backend=kernel), axis=1)
            transverse_blocks[block_index, q_index] = (
                to_numpy(trans_spec[freq_slice].real, kernel) / block_norm
            )

    return CurrentCorrelationResult(
        qpoints_cartesian=qpoints,
        frequencies_thz=freq_out,
        longitudinal=np.mean(longitudinal_blocks, axis=0),
        transverse=np.mean(transverse_blocks, axis=0),
        num_blocks=num_blocks,
        longitudinal_sem=standard_error(longitudinal_blocks),
        transverse_sem=standard_error(transverse_blocks),
        metadata={
            "dt_seconds": float(dt),
            "positive_only": bool(positive_only),
            "backend": kernel.name,
            "normalization": "block spectra divided by block_size * num_atoms, then averaged",
            "uncertainty": "standard error of the block-averaged spectra",
        },
    )

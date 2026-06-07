"""Kinematic q-EELS visibility and map construction."""

from dataclasses import dataclass

import numpy as np

from pySED.probe_weights import electron_form_factor_weights
from pySED.q_advisor import reciprocal_lattice, wrap_reduced


@dataclass
class EELSVisibility:
    qpoints_reduced: np.ndarray
    g_vectors_reduced: np.ndarray
    Q_reduced: np.ndarray
    Q_cartesian: np.ndarray
    visibility: np.ndarray
    metadata: dict


@dataclass
class EELSMap:
    Q_reduced: np.ndarray
    Q_cartesian: np.ndarray
    energy_axis: np.ndarray
    intensity: np.ndarray
    visibility: np.ndarray


@dataclass
class EELSDecompositionMap:
    Q_reduced: np.ndarray
    Q_cartesian: np.ndarray
    energy_axis: np.ndarray
    intensity: np.ndarray
    atom_intensity: np.ndarray
    direction_intensity: np.ndarray
    atom_interference_intensity: np.ndarray
    visibility: np.ndarray
    atom_visibility: np.ndarray
    direction_visibility: np.ndarray
    atom_interference: np.ndarray
    metadata: dict


@dataclass
class EELSVisibilityDecomposition:
    qpoints_reduced: np.ndarray
    g_vectors_reduced: np.ndarray
    Q_reduced: np.ndarray
    Q_cartesian: np.ndarray
    amplitude: np.ndarray
    visibility: np.ndarray
    atom_amplitude: np.ndarray
    atom_visibility: np.ndarray
    direction_amplitude: np.ndarray
    direction_visibility: np.ndarray
    atom_interference: np.ndarray
    metadata: dict


def lorentzian(axis, center, hwhm):
    axis = np.asarray(axis, dtype=float)
    return 1.0 / (1.0 + ((axis - center) / hwhm) ** 2)


def gaussian(axis, center, sigma):
    axis = np.asarray(axis, dtype=float)
    return np.exp(-0.5 * ((axis - center) / sigma) ** 2)


def fold_to_first_bz(qpoints_reduced):
    return wrap_reduced(qpoints_reduced, centered=True)


def _form_factors_for_qg(electron_form_factors, i_q, i_g, n_q, n_g, n_basis):
    if electron_form_factors is None:
        return np.ones(n_basis, dtype=float)

    form_arr = np.asarray(electron_form_factors, dtype=float)
    if form_arr.ndim == 1:
        if form_arr.shape[0] != n_basis:
            raise ValueError("1D electron_form_factors length must match number of basis atoms")
        return form_arr
    if form_arr.ndim == 2:
        if form_arr.shape == (n_q * n_g, n_basis):
            return form_arr[i_q * n_g + i_g]
        if form_arr.shape == (n_q, n_basis):
            return form_arr[i_q]
        raise ValueError(
            "2D electron_form_factors must have shape (num_q*num_g, num_basis) or (num_q, num_basis)"
        )
    if form_arr.ndim == 3:
        if form_arr.shape != (n_q, n_g, n_basis):
            raise ValueError("3D electron_form_factors must have shape (num_q, num_g, num_basis)")
        return form_arr[i_q, i_g]
    raise ValueError("electron_form_factors must be 1D, 2D, or 3D")


def _electron_form_factor_array(
    electron_form_factors,
    electron_form_factor_model,
    atom_types,
    qpoints_cartesian,
    coefficients_table,
    missing,
    n_q,
    n_g,
    n_basis,
):
    if electron_form_factors is not None:
        return electron_form_factors, "user supplied"

    model = "unit" if electron_form_factor_model is None else str(electron_form_factor_model).lower()
    model = model.replace("_", "-")
    if model in ("unit", "unity", "ones"):
        return None, "unit"
    if model not in ("mott-bethe", "mottbethe"):
        raise ValueError("electron_form_factor_model must be 'unit' or 'mott-bethe'")

    if atom_types is None:
        raise ValueError("atom_types are required for electron_form_factor_model='mott-bethe'")
    if len(atom_types) != n_basis:
        raise ValueError("atom_types length must match number of basis atoms")

    qpoints = np.asarray(qpoints_cartesian, dtype=float).reshape(n_q * n_g, 3)
    weights = electron_form_factor_weights(
        atom_types,
        qpoints,
        coefficients_table=coefficients_table,
        missing=missing,
    )
    return weights.reshape(n_q, n_g, n_basis), "mott-bethe"


def compute_mode_visibility(
    qpoints_reduced,
    g_vectors_reduced,
    primitive_cell,
    basis_positions_reduced,
    masses,
    eigenvectors,
    electron_form_factors=None,
    atom_types=None,
    electron_form_factor_model="unit",
    electron_form_factor_table=None,
    electron_missing="raise",
    pairwise_g_vectors=False,
):
    """Compute branch visibility for extended-zone kinematic q-EELS.

    The implemented factor is

        ``abs(sum_b f_b(Q) (Q dot e_bv(q)) exp(i 2*pi Q_red dot tau_b)
        / sqrt(m_b))**2``.
    """

    decomposition = compute_mode_visibility_decomposition(
        qpoints_reduced,
        g_vectors_reduced,
        primitive_cell,
        basis_positions_reduced,
        masses,
        eigenvectors,
        electron_form_factors=electron_form_factors,
        atom_types=atom_types,
        electron_form_factor_model=electron_form_factor_model,
        electron_form_factor_table=electron_form_factor_table,
        electron_missing=electron_missing,
        pairwise_g_vectors=pairwise_g_vectors,
    )
    return EELSVisibility(
        qpoints_reduced=decomposition.qpoints_reduced,
        g_vectors_reduced=decomposition.g_vectors_reduced,
        Q_reduced=decomposition.Q_reduced,
        Q_cartesian=decomposition.Q_cartesian,
        visibility=decomposition.visibility,
        metadata=decomposition.metadata,
    )


def compute_mode_visibility_decomposition(
    qpoints_reduced,
    g_vectors_reduced,
    primitive_cell,
    basis_positions_reduced,
    masses,
    eigenvectors,
    electron_form_factors=None,
    atom_types=None,
    electron_form_factor_model="unit",
    electron_form_factor_table=None,
    electron_missing="raise",
    pairwise_g_vectors=False,
):
    """Compute complex q-EELS visibility contributions.

    The primitive contribution is

    ``A[b, alpha] = f_b(Q) Q_alpha e_balpha(q, nu)
    exp(i 2*pi Q_red dot tau_b) / sqrt(m_b)``.

    Atom and direction amplitudes are coherent sums of ``A[b, alpha]``.  Their
    squared magnitudes are diagnostic quantities; they do not add up to the
    total visibility when there is destructive or constructive interference.
    """

    qpoints = np.atleast_2d(np.asarray(qpoints_reduced, dtype=float))
    g_vectors = np.atleast_2d(np.asarray(g_vectors_reduced, dtype=float))
    basis_pos = np.asarray(basis_positions_reduced, dtype=float)
    masses = np.asarray(masses, dtype=float)
    eig = np.asarray(eigenvectors, dtype=complex)

    if eig.shape[0] != qpoints.shape[0]:
        raise ValueError("eigenvectors first axis must match qpoints")
    if eig.shape[2:] != (basis_pos.shape[0], 3):
        raise ValueError("eigenvectors must have shape (num_q, num_modes, num_basis, 3)")
    if masses.shape[0] != basis_pos.shape[0]:
        raise ValueError("masses must contain one value per basis atom")

    reciprocal = reciprocal_lattice(primitive_cell)
    pairwise = bool(pairwise_g_vectors)
    if pairwise and g_vectors.shape != qpoints.shape:
        raise ValueError("pairwise g_vectors must have the same shape as qpoints")
    n_q = qpoints.shape[0]
    n_g = 1 if pairwise else g_vectors.shape[0]
    n_modes = eig.shape[1]
    n_basis = basis_pos.shape[0]
    Q_red = np.zeros((n_q, n_g, 3), dtype=float)
    Q_cart = np.zeros_like(Q_red)
    visibility = np.zeros((n_q, n_g, n_modes), dtype=float)
    amplitude = np.zeros((n_q, n_g, n_modes), dtype=complex)
    atom_amplitude = np.zeros((n_q, n_g, n_modes, n_basis), dtype=complex)
    direction_amplitude = np.zeros((n_q, n_g, n_modes, 3), dtype=complex)

    for i_q, q_red in enumerate(qpoints):
        active_g_vectors = g_vectors[i_q:i_q + 1] if pairwise else g_vectors
        for i_g, g_red in enumerate(active_g_vectors):
            big_q = q_red + g_red
            big_Q_cart = big_q @ reciprocal
            Q_red[i_q, i_g] = big_q
            Q_cart[i_q, i_g] = big_Q_cart

    electron_form_factor_data, electron_model = _electron_form_factor_array(
        electron_form_factors,
        electron_form_factor_model,
        atom_types,
        Q_cart,
        electron_form_factor_table,
        electron_missing,
        n_q,
        n_g,
        n_basis,
    )

    for i_q, q_red in enumerate(qpoints):
        active_g_vectors = g_vectors[i_q:i_q + 1] if pairwise else g_vectors
        for i_g, g_red in enumerate(active_g_vectors):
            big_q = q_red + g_red
            big_Q_cart = Q_cart[i_q, i_g]
            phase = np.exp(2j * np.pi * (basis_pos @ big_q))
            form = _form_factors_for_qg(
                electron_form_factor_data,
                i_q,
                i_g,
                n_q,
                n_g,
                n_basis,
            )

            prefactor = form * phase / np.sqrt(masses)
            component = eig[i_q] * big_Q_cart.reshape(1, 1, 3)
            contribution = component * prefactor.reshape(1, n_basis, 1)
            atom_amp = np.sum(contribution, axis=2)
            direction_amp = np.sum(contribution, axis=1)
            total_amp = np.sum(atom_amp, axis=1)
            atom_amplitude[i_q, i_g] = atom_amp
            direction_amplitude[i_q, i_g] = direction_amp
            amplitude[i_q, i_g] = total_amp
            visibility[i_q, i_g] = np.abs(total_amp) ** 2

    atom_visibility = np.abs(atom_amplitude) ** 2
    direction_visibility = np.abs(direction_amplitude) ** 2
    atom_interference = visibility - np.sum(atom_visibility, axis=3)

    return EELSVisibilityDecomposition(
        qpoints_reduced=qpoints,
        g_vectors_reduced=g_vectors,
        Q_reduced=Q_red,
        Q_cartesian=Q_cart,
        amplitude=amplitude,
        visibility=visibility,
        atom_amplitude=atom_amplitude,
        atom_visibility=atom_visibility,
        direction_amplitude=direction_amplitude,
        direction_visibility=direction_visibility,
        atom_interference=atom_interference,
        metadata={
            "model": "kinematic one-phonon visibility",
            "extended_zone": "Q = q + G",
            "electron_form_factor_model": electron_model,
            "pairwise_g_vectors": pairwise,
            "decomposition": "complex amplitudes; diagnostic intensities are not additive under interference",
        },
    )


def compute_mode_visibility_for_q_path(
    qpoints_reduced,
    g_vectors_reduced,
    primitive_cell,
    basis_positions_reduced,
    masses,
    eigenvectors,
    electron_form_factors=None,
    atom_types=None,
    electron_form_factor_model="unit",
    electron_form_factor_table=None,
    electron_missing="raise",
):
    """Compute EELS visibility for a pairwise experimental path.

    Unlike :func:`compute_mode_visibility`, this function interprets
    ``qpoints_reduced[i]`` and ``g_vectors_reduced[i]`` as one experimental
    point ``Q_i = q_i + G_i``.  The returned visibility has shape
    ``(num_points, 1, num_modes)`` so it remains compatible with
    :func:`build_eels_map_from_mode_spectra` without creating the full
    q-by-G Cartesian product.
    """

    return compute_mode_visibility(
        qpoints_reduced,
        g_vectors_reduced,
        primitive_cell,
        basis_positions_reduced,
        masses,
        eigenvectors,
        electron_form_factors=electron_form_factors,
        atom_types=atom_types,
        electron_form_factor_model=electron_form_factor_model,
        electron_form_factor_table=electron_form_factor_table,
        electron_missing=electron_missing,
        pairwise_g_vectors=True,
    )


def build_eels_map(
    visibility,
    mode_frequencies,
    energy_axis,
    broadening,
    broadening_kind="lorentzian",
):
    """Convert mode visibility into an energy-resolved q-EELS map."""

    mode_freq = np.asarray(mode_frequencies, dtype=float)
    energy = np.asarray(energy_axis, dtype=float)
    vis = visibility.visibility
    if mode_freq.shape != vis.shape[:1] + vis.shape[2:]:
        if mode_freq.shape == (vis.shape[0], vis.shape[2]):
            pass
        else:
            raise ValueError("mode_frequencies must have shape (num_q, num_modes)")

    intensity = np.zeros((vis.shape[0], vis.shape[1], energy.size), dtype=float)
    for i_q in range(vis.shape[0]):
        for i_g in range(vis.shape[1]):
            for mode in range(vis.shape[2]):
                center = mode_freq[i_q, mode]
                if broadening_kind == "lorentzian":
                    profile = lorentzian(energy, center, broadening)
                elif broadening_kind == "gaussian":
                    profile = gaussian(energy, center, broadening)
                else:
                    raise ValueError("broadening_kind must be 'lorentzian' or 'gaussian'")
                intensity[i_q, i_g] += vis[i_q, i_g, mode] * profile

    return EELSMap(
        Q_reduced=visibility.Q_reduced,
        Q_cartesian=visibility.Q_cartesian,
        energy_axis=energy,
        intensity=intensity,
        visibility=vis,
    )


def build_eels_map_from_mode_spectra(visibility, mode_spectra, energy_axis):
    """Weight branch-resolved spectral functions into an EELS intensity map.

    ``mode_spectra`` must have shape ``(num_energy, num_q, num_modes)``.  This
    matches :attr:`pySED.eigen_sed.EigenSEDResult.sed`, allowing branch-resolved
    MD spectra to be converted directly into extended-zone q-EELS maps:

    ``I(Q, energy) = sum_mode visibility(Q, mode) * A(q, mode, energy)``.
    """

    vis = np.asarray(visibility.visibility, dtype=float)
    spectra = np.asarray(mode_spectra, dtype=float)
    energy = np.asarray(energy_axis, dtype=float)

    if spectra.ndim != 3:
        raise ValueError("mode_spectra must have shape (num_energy, num_q, num_modes)")
    if spectra.shape[0] != energy.size:
        raise ValueError("energy_axis length must match mode_spectra first dimension")
    if spectra.shape[1:] != (vis.shape[0], vis.shape[2]):
        raise ValueError("mode_spectra must have shape (num_energy, num_q, num_modes)")

    intensity = np.einsum("qgm,eqm->qge", vis, spectra)

    return EELSMap(
        Q_reduced=visibility.Q_reduced,
        Q_cartesian=visibility.Q_cartesian,
        energy_axis=energy,
        intensity=intensity,
        visibility=vis,
    )


def build_eels_decomposition_map_from_mode_spectra(decomposition, mode_spectra, energy_axis):
    """Convert mode-resolved EELS diagnostics into energy-resolved maps.

    ``decomposition`` should come from
    :func:`compute_mode_visibility_decomposition`.  The returned maps preserve
    the non-additive diagnostics:

    ``total = sum_atom(atom_intensity) + atom_interference_intensity``.
    """

    vis = np.asarray(decomposition.visibility, dtype=float)
    atom_vis = np.asarray(decomposition.atom_visibility, dtype=float)
    direction_vis = np.asarray(decomposition.direction_visibility, dtype=float)
    atom_interference = np.asarray(decomposition.atom_interference, dtype=float)
    spectra = np.asarray(mode_spectra, dtype=float)
    energy = np.asarray(energy_axis, dtype=float)

    if spectra.ndim != 3:
        raise ValueError("mode_spectra must have shape (num_energy, num_q, num_modes)")
    if spectra.shape[0] != energy.size:
        raise ValueError("energy_axis length must match mode_spectra first dimension")
    if spectra.shape[1:] != (vis.shape[0], vis.shape[2]):
        raise ValueError("mode_spectra must have shape (num_energy, num_q, num_modes)")

    intensity = np.einsum("qgm,eqm->qge", vis, spectra)
    atom_intensity = np.einsum("qgmb,eqm->qgbe", atom_vis, spectra)
    direction_intensity = np.einsum("qgmd,eqm->qgde", direction_vis, spectra)
    atom_interference_intensity = np.einsum("qgm,eqm->qge", atom_interference, spectra)

    return EELSDecompositionMap(
        Q_reduced=decomposition.Q_reduced,
        Q_cartesian=decomposition.Q_cartesian,
        energy_axis=energy,
        intensity=intensity,
        atom_intensity=atom_intensity,
        direction_intensity=direction_intensity,
        atom_interference_intensity=atom_interference_intensity,
        visibility=vis,
        atom_visibility=atom_vis,
        direction_visibility=direction_vis,
        atom_interference=atom_interference,
        metadata={
            "source": "branch-resolved mode spectra",
            "decomposition": "diagnostic atom/direction maps; atom terms are not additive without interference",
            **decomposition.metadata,
        },
    )

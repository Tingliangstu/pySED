"""Mode-resolved one-phonon visibility for neutron and X-ray scattering."""

from dataclasses import dataclass

import numpy as np

from pySED.dsf import neutron_weights
from pySED.probe_weights import xray_form_factor_weights
from pySED.q_advisor import reciprocal_lattice


@dataclass
class OnePhononVisibility:
    """Coherent one-phonon mode visibility and diagnostic amplitudes."""

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
    scattering_weights: np.ndarray
    metadata: dict


@dataclass
class OnePhononMap:
    """Energy-resolved coherent one-phonon visibility map."""

    Q_reduced: np.ndarray
    Q_cartesian: np.ndarray
    energy_axis: np.ndarray
    intensity: np.ndarray
    visibility: np.ndarray
    metadata: dict


@dataclass
class OnePhononDecompositionMap:
    """Energy-resolved atom/direction/interference one-phonon diagnostics."""

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


def lorentzian(axis, center, hwhm):
    axis = np.asarray(axis, dtype=float)
    return 1.0 / (1.0 + ((axis - center) / hwhm) ** 2)


def gaussian(axis, center, sigma):
    axis = np.asarray(axis, dtype=float)
    return np.exp(-0.5 * ((axis - center) / sigma) ** 2)


def _coerce_weight_array(weights, n_q, n_g, n_basis):
    if weights is None:
        return np.ones((n_q, n_g, n_basis), dtype=float)
    arr = np.asarray(weights, dtype=float)
    if arr.ndim == 1:
        if arr.shape != (n_basis,):
            raise ValueError("1D scattering_weights length must match number of basis atoms")
        return np.tile(arr.reshape(1, 1, n_basis), (n_q, n_g, 1))
    if arr.ndim == 2:
        if arr.shape == (n_q * n_g, n_basis):
            return arr.reshape(n_q, n_g, n_basis)
        if arr.shape == (n_q, n_basis):
            return np.tile(arr[:, None, :], (1, n_g, 1))
        raise ValueError(
            "2D scattering_weights must have shape (num_q*num_g, num_basis) or (num_q, num_basis)"
        )
    if arr.ndim == 3:
        if arr.shape != (n_q, n_g, n_basis):
            raise ValueError("3D scattering_weights must have shape (num_q, num_g, num_basis)")
        return arr
    raise ValueError("scattering_weights must be 1D, 2D, or 3D")


def _frequency_scale(mode_frequencies, frequency_power, n_q, n_modes):
    if frequency_power is None or frequency_power == 0:
        return np.ones((n_q, n_modes), dtype=float)
    freq = np.asarray(mode_frequencies, dtype=float)
    if freq.shape != (n_q, n_modes):
        raise ValueError("mode_frequencies must have shape (num_q, num_modes)")
    scale = np.zeros_like(freq, dtype=float)
    positive = freq > 0.0
    scale[positive] = freq[positive] ** (0.5 * float(frequency_power))
    return scale


def _resolve_probe_weights(
    experiment,
    atom_types,
    qpoints_cartesian,
    scattering_weights,
    n_q,
    n_g,
    n_basis,
    xray_form_factor_table=None,
    xray_missing="atomic_number",
):
    if scattering_weights is not None:
        return _coerce_weight_array(scattering_weights, n_q, n_g, n_basis), "user supplied"

    experiment_key = str(experiment).lower()
    if experiment_key == "custom":
        return _coerce_weight_array(None, n_q, n_g, n_basis), "unit custom weights"
    if atom_types is None:
        raise ValueError("atom_types must be supplied for neutron or xray visibility")
    if len(atom_types) != n_basis:
        raise ValueError("atom_types length must match number of basis atoms")

    if experiment_key == "neutron":
        weights, _ = neutron_weights(atom_types, include_incoherent=False)
        return _coerce_weight_array(weights, n_q, n_g, n_basis), "neutron coherent scattering lengths"
    if experiment_key == "xray":
        flat_weights = xray_form_factor_weights(
            atom_types,
            qpoints_cartesian=qpoints_cartesian.reshape(-1, 3),
            coefficients_table=xray_form_factor_table,
            missing=xray_missing,
        )
        return _coerce_weight_array(flat_weights, n_q, n_g, n_basis), "Cromer-Mann X-ray form factors"

    raise ValueError("experiment must be 'neutron', 'xray', or 'custom'")


def compute_one_phonon_visibility(
    qpoints_reduced,
    g_vectors_reduced,
    primitive_cell,
    basis_positions_reduced,
    masses,
    eigenvectors,
    atom_types=None,
    experiment="neutron",
    scattering_weights=None,
    mode_frequencies=None,
    frequency_power=0,
    xray_form_factor_table=None,
    xray_missing="atomic_number",
    pairwise_g_vectors=False,
):
    """Compute coherent one-phonon INS/IXS visibility for each mode.

    The elementary complex contribution is

    ``A[b, alpha] = w_b(Q) Q_alpha e_balpha(q, nu)
    exp(i 2*pi Q_red dot tau_b) / sqrt(m_b)``.

    ``frequency_power=-1`` adds the common harmonic one-phonon ``1/omega``
    intensity prefactor.  The default ``0`` reports the pure polarization,
    probe-weight, mass, and interference visibility.
    """

    qpoints = np.atleast_2d(np.asarray(qpoints_reduced, dtype=float))
    g_vectors = np.atleast_2d(np.asarray(g_vectors_reduced, dtype=float))
    basis_pos = np.asarray(basis_positions_reduced, dtype=float)
    masses = np.asarray(masses, dtype=float)
    eig = np.asarray(eigenvectors, dtype=complex)

    if qpoints.shape[1] != 3 or g_vectors.shape[1] != 3:
        raise ValueError("qpoints_reduced and g_vectors_reduced must have shape (n, 3)")
    if basis_pos.ndim != 2 or basis_pos.shape[1] != 3:
        raise ValueError("basis_positions_reduced must have shape (num_basis, 3)")
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
    for i_q, q_red in enumerate(qpoints):
        active_g_vectors = g_vectors[i_q:i_q + 1] if pairwise else g_vectors
        for i_g, g_red in enumerate(active_g_vectors):
            Q_red[i_q, i_g] = q_red + g_red
    Q_cart = Q_red @ reciprocal

    weights, weight_model = _resolve_probe_weights(
        experiment,
        atom_types,
        Q_cart,
        scattering_weights,
        n_q,
        n_g,
        n_basis,
        xray_form_factor_table=xray_form_factor_table,
        xray_missing=xray_missing,
    )
    freq_scale = _frequency_scale(mode_frequencies, frequency_power, n_q, n_modes)

    amplitude = np.zeros((n_q, n_g, n_modes), dtype=complex)
    visibility = np.zeros((n_q, n_g, n_modes), dtype=float)
    atom_amplitude = np.zeros((n_q, n_g, n_modes, n_basis), dtype=complex)
    direction_amplitude = np.zeros((n_q, n_g, n_modes, 3), dtype=complex)

    for i_q in range(n_q):
        for i_g in range(n_g):
            phase = np.exp(2j * np.pi * (basis_pos @ Q_red[i_q, i_g]))
            prefactor = weights[i_q, i_g] * phase / np.sqrt(masses)
            component = eig[i_q] * Q_cart[i_q, i_g].reshape(1, 1, 3)
            contribution = component * prefactor.reshape(1, n_basis, 1)
            contribution *= freq_scale[i_q].reshape(n_modes, 1, 1)
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

    return OnePhononVisibility(
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
        scattering_weights=weights,
        metadata={
            "model": "coherent one-phonon neutron/xray visibility",
            "experiment": str(experiment).lower(),
            "weight_model": weight_model,
            "extended_zone": "Q = q + G",
            "frequency_power": frequency_power,
            "pairwise_g_vectors": pairwise,
            "decomposition": "complex amplitudes; diagnostic intensities are not additive under interference",
        },
    )


def compute_one_phonon_visibility_for_q_path(
    qpoints_reduced,
    g_vectors_reduced,
    primitive_cell,
    basis_positions_reduced,
    masses,
    eigenvectors,
    atom_types=None,
    experiment="neutron",
    scattering_weights=None,
    mode_frequencies=None,
    frequency_power=0,
    xray_form_factor_table=None,
    xray_missing="atomic_number",
):
    """Compute coherent INS/IXS one-phonon visibility for a pairwise Q path."""

    return compute_one_phonon_visibility(
        qpoints_reduced,
        g_vectors_reduced,
        primitive_cell,
        basis_positions_reduced,
        masses,
        eigenvectors,
        atom_types=atom_types,
        experiment=experiment,
        scattering_weights=scattering_weights,
        mode_frequencies=mode_frequencies,
        frequency_power=frequency_power,
        xray_form_factor_table=xray_form_factor_table,
        xray_missing=xray_missing,
        pairwise_g_vectors=True,
    )


def build_one_phonon_map(
    visibility,
    mode_frequencies,
    energy_axis,
    broadening,
    broadening_kind="lorentzian",
):
    """Convert one-phonon visibility into a broadened harmonic map."""

    vis = np.asarray(visibility.visibility, dtype=float)
    mode_freq = np.asarray(mode_frequencies, dtype=float)
    energy = np.asarray(energy_axis, dtype=float)
    if mode_freq.shape != (vis.shape[0], vis.shape[2]):
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

    return OnePhononMap(
        Q_reduced=visibility.Q_reduced,
        Q_cartesian=visibility.Q_cartesian,
        energy_axis=energy,
        intensity=intensity,
        visibility=vis,
        metadata={
            "source": "harmonic mode frequencies",
            "broadening": float(broadening),
            "broadening_kind": broadening_kind,
            **visibility.metadata,
        },
    )


def build_one_phonon_map_from_mode_spectra(visibility, mode_spectra, energy_axis):
    """Weight branch-resolved spectra into coherent INS/IXS intensity maps.

    ``mode_spectra`` must have shape ``(num_energy, num_q, num_modes)`` and can
    be taken directly from :attr:`pySED.eigen_sed.EigenSEDResult.sed`.
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
    return OnePhononMap(
        Q_reduced=visibility.Q_reduced,
        Q_cartesian=visibility.Q_cartesian,
        energy_axis=energy,
        intensity=intensity,
        visibility=vis,
        metadata={
            "source": "branch-resolved mode spectra",
            **visibility.metadata,
        },
    )


def build_one_phonon_decomposition_map_from_mode_spectra(visibility, mode_spectra, energy_axis):
    """Build energy-resolved atom/direction/interference INS/IXS maps."""

    vis = np.asarray(visibility.visibility, dtype=float)
    atom_vis = np.asarray(visibility.atom_visibility, dtype=float)
    direction_vis = np.asarray(visibility.direction_visibility, dtype=float)
    atom_interference = np.asarray(visibility.atom_interference, dtype=float)
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

    return OnePhononDecompositionMap(
        Q_reduced=visibility.Q_reduced,
        Q_cartesian=visibility.Q_cartesian,
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
            **visibility.metadata,
        },
    )

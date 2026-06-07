"""Eigenvector-resolved spectral energy density."""

from dataclasses import dataclass

import numpy as np

from pySED.q_advisor import QAdvisor
from pySED.scattering_kernel import frequency_grid, resolve_backend, to_numpy


@dataclass
class EigenvectorSet:
    qpoints_reduced: np.ndarray
    frequencies: np.ndarray
    eigenvectors: np.ndarray
    source: str = "unknown"


@dataclass
class EigenSEDResult:
    frequencies_thz: np.ndarray
    qpoints_reduced: np.ndarray
    sed: np.ndarray
    mode_frequencies: np.ndarray
    metadata: dict


def read_phonopy_band_yaml(path):
    """Read q points, frequencies, and eigenvectors from a phonopy band.yaml."""

    try:
        import yaml
    except ImportError as exc:
        raise ImportError("PyYAML is required to read phonopy band.yaml files") from exc

    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle)

    phonons = data.get("phonon", [])
    qpoints = []
    frequencies = []
    eigenvectors = []
    for phonon in phonons:
        qpoints.append(phonon["q-position"])
        mode_freqs = []
        mode_eigs = []
        for band in phonon["band"]:
            mode_freqs.append(band.get("frequency", np.nan))
            eig = []
            for atom in band["eigenvector"]:
                atom_vec = []
                for component in atom:
                    atom_vec.append(complex(component[0], component[1]))
                eig.append(atom_vec)
            mode_eigs.append(eig)
        frequencies.append(mode_freqs)
        eigenvectors.append(mode_eigs)

    return EigenvectorSet(
        qpoints_reduced=np.asarray(qpoints, dtype=float),
        frequencies=np.asarray(frequencies, dtype=float),
        eigenvectors=np.asarray(eigenvectors, dtype=complex),
        source=str(path),
    )


def validate_eigenvector_qpoints(eigenvectors, advisor, target_qpoints=None, atol=1e-6):
    qpoints = np.asarray(eigenvectors.qpoints_reduced, dtype=float)
    if not advisor.validate_qpoints(qpoints, atol=atol):
        raise ValueError("eigenvector q points are not commensurate with the MD supercell")

    if target_qpoints is not None:
        target = np.asarray(target_qpoints, dtype=float)
        if target.shape != qpoints.shape:
            raise ValueError("target q points and eigenvector q points have different shapes")
        delta = qpoints - target
        delta -= np.round(delta)
        if not np.allclose(delta, 0.0, atol=atol):
            raise ValueError("eigenvector q points do not match the requested pySED q points")
    return True


def _basis_atom_ids_from_index(basis_index):
    basis_index = np.asarray(basis_index, dtype=int)
    labels = np.unique(basis_index)
    if not np.array_equal(labels, np.arange(1, labels.size + 1)):
        raise ValueError("basis_index must use contiguous one-based labels")
    groups = [np.flatnonzero(basis_index == label) for label in labels]
    lengths = {len(group) for group in groups}
    if len(lengths) != 1:
        raise ValueError("each basis atom must appear once in every unit cell")
    return np.asarray(groups, dtype=int)


def compute_eigen_sed(
    velocities,
    unitcell_vectors,
    basis_index,
    masses,
    primitive_cell,
    supercell_cell,
    eigenvectors,
    dt,
    num_blocks=1,
    target_qpoints=None,
    validate_q=True,
    backend="cpu",
):
    """Project velocities onto phonon eigenvectors and compute branch SED.

    ``velocities`` must have shape ``(num_frames, num_atoms, 3)``.  The
    ``basis_index`` array maps each atom to a one-based primitive basis label.
    ``unitcell_vectors`` must contain the equilibrium cell-vector position of
    each repeated primitive cell in Angstrom.
    """

    velocities = np.asarray(velocities, dtype=float)
    unitcell_vectors = np.asarray(unitcell_vectors, dtype=float)
    basis_ids = _basis_atom_ids_from_index(basis_index)
    masses = np.asarray(masses, dtype=float)
    kernel = resolve_backend(backend)
    xp = kernel.xp

    if velocities.ndim != 3 or velocities.shape[2] != 3:
        raise ValueError("velocities must have shape (num_frames, num_atoms, 3)")
    if masses.shape[0] != basis_ids.shape[0]:
        raise ValueError("masses must contain one value per basis atom")
    if unitcell_vectors.shape != (basis_ids.shape[1], 3):
        raise ValueError("unitcell_vectors must have shape (num_unit_cells, 3)")

    advisor = QAdvisor(primitive_cell, supercell_cell)
    if validate_q:
        validate_eigenvector_qpoints(eigenvectors, advisor, target_qpoints=target_qpoints)

    qpoints_red = np.asarray(eigenvectors.qpoints_reduced, dtype=float)
    qpoints_cart = advisor.reduced_to_cartesian(qpoints_red)
    eig = np.asarray(eigenvectors.eigenvectors, dtype=complex)
    if eig.shape[:1] != (qpoints_red.shape[0],):
        raise ValueError("eigenvectors first axis must match qpoints")
    if eig.shape[2:] != (basis_ids.shape[0], 3):
        raise ValueError("eigenvectors must have shape (num_q, num_modes, num_basis, 3)")

    n_frames = velocities.shape[0]
    block_size = n_frames // num_blocks
    if block_size < 2:
        raise ValueError("each block must contain at least two frames")
    grid = frequency_grid(block_size, dt, positive_only=True, backend=kernel)
    pos_slice = grid.freq_slice
    freq_out = grid.frequencies_thz

    n_q, n_modes = eig.shape[0], eig.shape[1]
    sed = xp.zeros((len(freq_out), n_q, n_modes), dtype=float)

    velocities_backend = xp.asarray(velocities, dtype=float)
    unitcell_vectors_backend = xp.asarray(unitcell_vectors, dtype=float)
    qpoints_cart_backend = xp.asarray(qpoints_cart, dtype=float)
    eig_backend = xp.asarray(eig, dtype=complex)
    sqrt_m = xp.sqrt(xp.asarray(masses, dtype=float))
    for block in range(num_blocks):
        start = block * block_size
        stop = start + block_size
        vel_block = velocities_backend[start:stop]

        for q_index, q_cart in enumerate(qpoints_cart):
            q_cart_backend = qpoints_cart_backend[q_index]
            phase = xp.exp(1j * (unitcell_vectors_backend @ q_cart_backend))
            basis_fft_inputs = xp.zeros((block_size, basis_ids.shape[0], 3), dtype=complex)
            for b_idx, atom_ids in enumerate(basis_ids):
                basis_vel = vel_block[:, atom_ids, :]
                phased = basis_vel * phase.reshape(1, -1, 1)
                basis_fft_inputs[:, b_idx, :] = xp.sum(phased, axis=1) * sqrt_m[b_idx]

            projected = xp.einsum(
                "tba,mba->tm",
                basis_fft_inputs,
                xp.conjugate(eig_backend[q_index]),
            )
            spectrum = xp.abs(kernel.fft(projected, axis=0)) ** 2
            sed[:, q_index, :] += xp.real(spectrum[pos_slice])

    sed = to_numpy(sed / float(num_blocks), kernel)
    return EigenSEDResult(
        frequencies_thz=freq_out,
        qpoints_reduced=qpoints_red,
        sed=sed,
        mode_frequencies=np.asarray(eigenvectors.frequencies, dtype=float),
        metadata={
            "dt_seconds": float(dt),
            "num_blocks": int(num_blocks),
            "backend": kernel.name,
            "normalization": "branch projected |FFT(qdot)|^2 averaged over blocks",
        },
    )

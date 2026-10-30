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


"""SED mathematics: spatial projection, time FFT, and normalization."""

from dataclasses import dataclass

import numpy as np
from scipy.fftpack import fft


JOULE_PER_EV = 1.602176634e-19
JOULE_SECOND_TO_EV_PER_THZ = 2 * np.pi * 1.0e12 / JOULE_PER_EV


@dataclass(frozen=True)
class SEDKernelConfig:
    """Static topology shared by the mathematical kernel and both backends."""

    basis_ids: tuple
    masses: np.ndarray
    basis_to_type: tuple
    qpoints: np.ndarray
    num_types: int
    output_partial: bool
    amu_2_kg: float


def calculate_q_batch(
    vels,
    cell_vecs,
    qpoints,
    basis_ids_by_basis,
    config,
    xp=np,
    fft_function=fft,
    project_function=None,
):
    """Compute unnormalized SED; this function owns no files or workers.

    Shapes: vels (T, atoms, 3), cell_vecs (cells, 3), qpoints (Q, 3).
    Each basis selects one atom per cell. The projection sums only cells,
    producing (T, 3, Q); FFT acts only on T. The result is (T, Q), or
    (T, Q, types, 3) for partial SED. Here T becomes the full FFT grid.

    NumPy and CuPy use the same operations. Backends may wrap projection
    and FFT for timing; the physics does not depend on those wrappers.
    """
    project = xp.tensordot if project_function is None else project_function

    # 1. Spatial phase for every cell and Q point: (cells, Q).
    phases = xp.exp(1.0j * (cell_vecs @ qpoints.T))
    num_frames = int(vels.shape[0])
    num_qpoints = int(qpoints.shape[0])
    shape = (num_frames, num_qpoints)
    if config.output_partial:
        shape += (config.num_types, 3)
    result = xp.zeros(shape, dtype=xp.float64)

    for b_idx, basis_ids in enumerate(basis_ids_by_basis):
        # 2. Sum velocities over cells BEFORE FFT (FFT is linear).
        projected_vels = project(
            vels[:, basis_ids, :], phases, axes=([1], [0]),
        )
        # 3. Time FFT, then mass-weighted power: (frequency, xyz, Q).
        transformed_vels = fft_function(projected_vels, axis=0)
        spectra = (
            xp.abs(transformed_vels) ** 2
            * config.masses[b_idx]
            * config.amu_2_kg
        )
        del projected_vels, transformed_vels

        # 4. Retain element/xyz components or sum them into the total SED.
        if config.output_partial:
            type_index = config.basis_to_type[b_idx]
            result[:, :, type_index, :] += spectra.transpose(0, 2, 1)
        else:
            result += spectra.sum(axis=1)
    return result


def block_normalization(dt, duration, num_unit_cells):
    """Scale FFT power; dt and duration are seconds, velocities are m/s."""
    return dt**2 / (4 * np.pi * duration * num_unit_cells)


def finalize_spectrum(sed_sum, frequencies, num_blocks):
    """Average and fold in place; return one-sided SED (eV/THz) and THz.

    Keep the existing 3.0 frequency selection, including exclusion of the
    Nyquist bin. Do not change this convention as part of a code move.
    """
    sed_sum /= num_blocks
    n_half = len(frequencies) // 2
    negative = sed_sum[:n_half:-1, ...]
    paired_length = negative.shape[0]
    positive = sed_sum[1 : 1 + paired_length, ...]
    sed_sum[1 : 1 + paired_length, ...] = positive + negative
    sed_avg = sed_sum[:n_half, ...]
    sed_avg *= JOULE_SECOND_TO_EV_PER_THZ
    return sed_avg, frequencies[:n_half]

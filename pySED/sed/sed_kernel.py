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
# The angular-frequency density gains 2*pi when expressed per ordinary
# frequency; 1e12 converts per Hz to per THz, and JOULE_PER_EV converts energy.
JOULE_SECOND_TO_EV_PER_THZ = 2 * np.pi * 1.0e12 / JOULE_PER_EV
SED_UNIT = "eV/THz"
LEGACY_SED_UNIT = "legacy 2.x scale"
LEGACY_UNIT_NOTICE = (
    "legacy_unit is for comparison with pySED 2.x results and will be "
    "removed in a future major release. Use legacy_unit = False for new calculations."
)


@dataclass(frozen=True)
class SEDKernelConfig:
    """Static atom mapping and q points shared by CPU and GPU calculations.

    ``basis_ids`` contains one zero-based atom-index array per basis atom,
    ordered consistently across cells. ``masses`` has one mass in amu per
    basis atom, and ``basis_to_type`` selects its partial-output channel.
    ``qpoints`` has shape (Q, 3) in inverse Angstrom with the 2*pi factor
    included; ``amu_2_kg`` converts masses for the kinetic-energy weighting.
    """

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
    """Compute mass-weighted FFT power for a batch of q points.

    Shapes: vels (T, atoms, 3), cell_vecs (cells, 3), qpoints (Q, 3).
    Velocities are in m/s; cell vectors are reference positions in Angstrom.
    Cartesian q points include 2*pi, making ``cell_vecs @ qpoints.T`` a phase
    in radians. ``basis_ids_by_basis`` selects matching atoms in each cell.
    Each basis selects one atom per cell. The projection sums only cells,
    producing (T, 3, Q); FFT acts only on T. The result is (T, Q), or
    (T, Q, types, 3) for partial SED. Here T becomes the full FFT grid.

    The same operations support NumPy and CuPy arrays through ``xp`` and
    ``fft_function``. This returns a new array of power; sampling-time
    normalization, block averaging, and one-sided folding are separate steps.
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


def block_normalization(dt, duration, num_unit_cells, *, legacy_unit=False):
    """Return the factor applied to each block's mass-weighted FFT power.

    ``dt`` is the saved-frame interval and ``duration`` is the block duration,
    both in seconds; ``num_unit_cells`` is the number of repeated cells.
    The physical factor is dt**2 / (4*pi*duration*num_unit_cells): the squared
    interval converts the squared discrete FFT sum into a squared time integral.
    Legacy 2.x scaling omits dt**2. Output-unit conversion is done later.
    """
    time_factor = 1.0 if legacy_unit else dt**2
    return time_factor / (4 * np.pi * duration * num_unit_cells)


def finalize_spectrum(sed_sum, frequencies, num_blocks, *, legacy_unit=False):
    """Average blocks and fold to eV/THz or the historical 2.x output scale.

    DC and all positive frequencies are retained; an even-length Nyquist
    bin is excluded. The default scale sums each positive/negative pair and
    converts to eV/THz. Legacy scaling averages each pair and skips that unit
    conversion. Frequencies remain in THz in both modes.

    ``sed_sum`` already contains block-normalized spectra, with frequency
    on axis 0 and optional q/type/direction axes afterward. ``frequencies``
    is the corresponding full FFT grid in THz. The input array is modified
    in place; the returned (spectrum, frequency_grid) uses slices of the
    original arrays. Use the same ``legacy_unit`` value as block_normalization.
    """
    sed_sum /= num_blocks
    n_half = len(frequencies) // 2
    n_positive = (len(frequencies) + 1) // 2
    negative = sed_sum[:n_half:-1, ...]
    paired_length = negative.shape[0]
    positive = sed_sum[1 : 1 + paired_length, ...]
    sed_sum[1 : 1 + paired_length, ...] = positive + negative
    sed_avg = sed_sum[:n_positive, ...]
    if legacy_unit:
        sed_avg[1:, ...] /= 2.0
    else:
        sed_avg *= JOULE_SECOND_TO_EV_PER_THZ
    return sed_avg, frequencies[:n_positive]

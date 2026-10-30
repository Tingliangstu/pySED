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


"""Backend-only helpers: Q batches, progress, and memory estimates."""

import math

import numpy as np


def qpoint_batches(num_qpoints, num_workers):
    """Split typical 20-40 Q points into small, balanced worker batches."""

    if num_qpoints < 1:
        return []
    if num_qpoints == 1:
        return [(0, 1)]

    # Choose enough batches to cap each at five Q points, but not so many that
    # a batch contains only one. Evenly distributing the remainder avoids a
    # short final batch.
    minimum_batches = math.ceil(num_qpoints / 5)
    maximum_batches = num_qpoints // 2
    num_batches = max(
        minimum_batches,
        min(num_workers, maximum_batches),
    )
    base_size, extra = divmod(num_qpoints, num_batches)

    batches = []
    start = 0
    for batch_index in range(num_batches):
        batch_size = base_size + (batch_index < extra)
        stop = start + batch_size
        batches.append((start, stop))
        start = stop
    return batches


def estimate_kernel_working_bytes(vels, config, max_q_batch):
    """Estimate the largest temporary arrays used by one SED kernel."""

    num_frames = vels.shape[0]
    max_basis_atoms = max(len(ids) for ids in config.basis_ids)
    basis_velocities = (
        num_frames * max_basis_atoms * 3 * vels.dtype.itemsize
    )
    phases = max_basis_atoms * max_q_batch * np.dtype(np.complex128).itemsize
    spectral_values = num_frames * 3 * max_q_batch
    transforms = spectral_values * (
        2 * np.dtype(np.complex128).itemsize
        + np.dtype(np.float64).itemsize
    )
    result_components = config.num_types * 3 if config.output_partial else 1
    result = (
        num_frames
        * max_q_batch
        * result_components
        * np.dtype(np.float64).itemsize
    )
    return basis_velocities + phases + transforms + result


def print_qpoints(reduced_qpoints, q_start, q_stop):
    """Keep the familiar one-line progress message for each Q point."""
    num_qpoints = len(reduced_qpoints)
    for q_index in range(q_start, q_stop):
        qpoint = reduced_qpoints[q_index]
        print(
            "\tNow calculating q-point "
            f"{q_index + 1}/{num_qpoints}:\t"
            f"q = ({qpoint[0]:.4f}, {qpoint[1]:.4f}, "
            f"{qpoint[2]:.4f})"
        )


def print_cpu_memory_estimate(
    vels, cell_vecs, config, accumulator, result, *,
    backend, num_workers=1, parallel=False, prefetch=False,
):
    """Report host-array memory without changing any calculation arrays."""
    batches = qpoint_batches(config.qpoints.shape[0], num_workers)
    max_q_batch = max(stop - start for start, stop in batches)
    worker_peak = estimate_kernel_working_bytes(vels, config, max_q_batch)
    persistent_arrays = (
        accumulator.nbytes + result.nbytes + vels.nbytes + cell_vecs.nbytes
    )
    prefetched_block = 2 * vels.nbytes + cell_vecs.nbytes if prefetch else 0
    shared_velocities = vels.nbytes if parallel else 0
    active_workers = 0 if backend == "cupy" else (num_workers if parallel else 1)
    estimated_bytes = (
        persistent_arrays + prefetched_block + shared_velocities
        + active_workers * worker_peak
    )
    print(
        "\n  Estimated peak CPU array memory: "
        f"{estimated_bytes / 1e9:.2f} GB "
        "(trajectory block reading + SED arrays; excluding Python overhead)."
    )

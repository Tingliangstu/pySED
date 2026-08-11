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

import math
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import contextmanager
from dataclasses import dataclass
from itertools import chain
from multiprocessing import shared_memory

import numpy as np
from scipy.fftpack import fft, fftfreq

from pySED.io.schema import POSITIONS, VELOCITIES
from pySED.structure.atoms import atomic_masses


JOULE_PER_EV = 1.602176634e-19
JOULE_SECOND_TO_EV_PER_THZ = (
    2 * np.pi * 1.0e12 / JOULE_PER_EV
)


@dataclass(frozen=True)
class _SEDKernelConfig:
    """Store the small, static arrays needed by every SED worker."""

    basis_ids: tuple
    masses: np.ndarray
    basis_to_type: tuple
    qpoints: np.ndarray
    num_types: int
    output_partial: bool
    amu_2_kg: float


class _GPUEventProfiler:
    """Accumulate asynchronous CUDA event timings with minimal interference."""

    def __init__(self, cupy_module):
        self.cp = cupy_module
        self.events = {}

    def start(self):
        """Record the beginning of one GPU stage."""

        event = self.cp.cuda.Event()
        event.record()
        return event

    def stop(self, stage, start_event):
        """Record the end of a stage without synchronizing the GPU."""

        end_event = self.cp.cuda.Event()
        end_event.record()
        self.events.setdefault(stage, []).append(
            (start_event, end_event)
        )

    def measure(self, stage, operation):
        """Enclose one queued GPU operation between CUDA events."""

        start_event = self.start()
        result = operation()
        self.stop(stage, start_event)
        return result

    def elapsed_seconds(self):
        """Synchronize once and return accumulated CUDA time by stage."""

        self.cp.cuda.get_current_stream().synchronize()
        return {
            stage: sum(
                self.cp.cuda.get_elapsed_time(start, end)
                for start, end in event_pairs
            )
            / 1000.0
            for stage, event_pairs in self.events.items()
        }


_CPU_WORKER_CONFIG = None


def _initialize_cpu_worker(config):
    """Install static topology metadata once in each persistent CPU worker."""

    global _CPU_WORKER_CONFIG
    _CPU_WORKER_CONFIG = config


def _calculate_q_batch(
    vels,
    cell_vecs,
    qpoints,
    basis_ids_by_basis,
    config,
    xp=np,
    fft_function=fft,
    gpu_profiler=None,
):
    """Calculate a Q-point batch with either NumPy or CuPy operations."""

    # One matrix multiplication creates phase factors for the whole Q batch.
    phases = xp.exp(1.0j * (cell_vecs @ qpoints.T))
    num_frames = int(vels.shape[0])
    num_qpoints = int(qpoints.shape[0])

    if config.output_partial:
        result = xp.zeros(
            (num_frames, num_qpoints, config.num_types, 3),
            dtype=xp.float64,
        )
    else:
        result = xp.zeros((num_frames, num_qpoints), dtype=xp.float64)

    for b_idx, basis_ids in enumerate(basis_ids_by_basis):
        # Contract only the unit-cell axis. Time, xyz, and Q remain separate.
        if gpu_profiler is None:
            projected_vels = xp.tensordot(
                vels[:, basis_ids, :],
                phases,
                axes=([1], [0]),
            )
        else:
            projected_vels = gpu_profiler.measure(
                "projection",
                lambda: xp.tensordot(
                    vels[:, basis_ids, :],
                    phases,
                    axes=([1], [0]),
                ),
            )

        if gpu_profiler is None:
            transformed_vels = fft_function(projected_vels, axis=0)
        else:
            transformed_vels = gpu_profiler.measure(
                "fft",
                lambda values=projected_vels: fft_function(values, axis=0),
            )
        spectra = (
            xp.abs(transformed_vels) ** 2
            * config.masses[b_idx]
            * config.amu_2_kg
        )
        del projected_vels, transformed_vels

        if config.output_partial:
            type_index = config.basis_to_type[b_idx]
            result[:, :, type_index, :] += spectra.transpose(0, 2, 1)
        else:
            result += spectra.sum(axis=1)

    return result


def _compute_shared_q_batch(
    shared_name,
    shape,
    dtype_string,
    cell_vecs,
    q_start,
    q_stop,
):
    """Attach to one shared velocity block and calculate one Q range."""

    if _CPU_WORKER_CONFIG is None:
        raise RuntimeError("SED worker was not initialized")

    segment = shared_memory.SharedMemory(name=shared_name)
    try:
        # The worker creates only an ndarray view; it does not copy vels.
        vels = np.ndarray(
            shape,
            dtype=np.dtype(dtype_string),
            buffer=segment.buf,
        )
        vels.setflags(write=False)
        result = _calculate_q_batch(
            vels=vels,
            cell_vecs=cell_vecs,
            qpoints=_CPU_WORKER_CONFIG.qpoints[q_start:q_stop],
            basis_ids_by_basis=_CPU_WORKER_CONFIG.basis_ids,
            config=_CPU_WORKER_CONFIG,
        )
        return q_start, q_stop, result
    finally:
        segment.close()


@contextmanager
def _share_numpy_array(values):
    """Copy one contiguous array into shared memory and always clean it up."""

    contiguous = np.ascontiguousarray(values)
    segment = shared_memory.SharedMemory(create=True, size=contiguous.nbytes)
    shared_values = np.ndarray(
        contiguous.shape,
        dtype=contiguous.dtype,
        buffer=segment.buf,
    )
    np.copyto(shared_values, contiguous)

    try:
        yield segment.name, contiguous.shape, contiguous.dtype.str
    finally:
        segment.close()
        segment.unlink()


def _qpoint_batches(num_qpoints, num_workers):
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


def _load_cupy():
    """Import CuPy lazily and verify that a CUDA device is usable."""

    try:
        import cupy as cp
    except ImportError as error:
        raise ImportError(
            "backend = cupy requires a CuPy package matching the installed "
            "CUDA version"
        ) from error

    try:
        device_count = cp.cuda.runtime.getDeviceCount()
    except Exception as error:
        raise RuntimeError(
            "CuPy is installed, but no usable CUDA device was found"
        ) from error
    if device_count < 1:
        raise RuntimeError("CuPy found no usable CUDA device")
    return cp


def _element_symbol_from_mass(mass, tolerance=5.0e-2):
    """Infer an element symbol from a standard mass without guessing."""

    candidates = [
        (abs(mass - standard_mass), symbol)
        for symbol, standard_mass in atomic_masses.items()
        if standard_mass is not None
        and abs(mass - standard_mass) <= tolerance
    ]
    if not candidates:
        return "unknown"
    return min(candidates)[1]


def _estimate_kernel_working_bytes(vels, config, max_q_batch):
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


class spectral_energy_density:
    """Calculate spectral energy density with a CPU or CuPy backend."""

    def __init__(self, params, trajectory_info):
        """Validate the trajectory and prepare the common frequency grid."""

        self.num_frame = (
            params.total_num_steps // params.output_data_stride
        )
        self.num_frames_per_block = self.num_frame // params.num_blocks

        if trajectory_info.n_atoms != params.num_atoms:
            raise ValueError(
                f"Trajectory has {trajectory_info.n_atoms} atoms, "
                f"but num_atoms = {params.num_atoms}."
            )

        print(
            "\nThe number of frames used to calculate SED is "
            f"{self.num_frame}."
        )
        print(
            f"\nThe {self.num_frame} frame trajectories were divided into "
            f"{params.num_blocks} blocks for averaging."
        )

        # Convert the saved-frame interval from fs to seconds.
        self.dt = params.time_step * params.output_data_stride / 1e15

        # Total simulated time represented by one averaging block.
        self.t_o = (
            params.time_step
            * params.total_num_steps
            / params.num_blocks
            / 1e15
        )
        self.amu_2_kg = 1.66054e-27

        print(
            "\n****** Velocity unit normalized to m/s; "
            "SED output unit is eV/THz ******"
        )

        # scipy.fftpack.fftfreq returns Hz; expose the frequency axis in THz.
        self.freq_fft = (
            fftfreq(self.num_frames_per_block, self.dt) / 1e12
        )

    def compute_sed(
        self,
        params,
        lattice_info,
        first_trajectory_block,
        remaining_trajectory_blocks,
    ):
        """Prepare static metadata and run every block through one backend."""

        self.backend = params.backend
        self.trajectory_prefetch = bool(params.trajectory_prefetch)
        max_cores = params.max_cores
        start_time = time.time()

        # These topology values do not change between blocks or Q points.
        self.num_unit_cells = lattice_info.unitcell_index.max()
        num_basis = lattice_info.basis_index.max()
        self.num_blocks = params.num_blocks

        unique_masses = np.unique(lattice_info.masses)
        self.num_types = len(unique_masses)
        type_symbols = []
        for type_index, mass in enumerate(unique_masses):
            symbol = _element_symbol_from_mass(mass)
            if symbol == "unknown":
                symbol = f"type{type_index + 1}"
            type_symbols.append(symbol)
        self.type_symbols = tuple(type_symbols)
        basis_masses = np.asarray(
            lattice_info.masses[:num_basis],
            dtype=float,
        )
        basis_to_type = [
            np.where(unique_masses == mass)[0][0]
            for mass in basis_masses
        ]
        self.output_partial = getattr(params, "output_partial", 0)

        # Basis membership is topology-only, so search the atoms only once.
        basis_ids = tuple(
            np.flatnonzero(lattice_info.basis_index == (b_idx + 1))
            for b_idx in range(num_basis)
        )
        for b_idx, atom_ids in enumerate(basis_ids):
            if atom_ids.size != self.num_unit_cells:
                raise ValueError(
                    f"Basis {b_idx + 1} contains {atom_ids.size} atoms; "
                    f"expected {self.num_unit_cells} unit cells."
                )

        kernel_config = _SEDKernelConfig(
            basis_ids=basis_ids,
            masses=basis_masses,
            basis_to_type=tuple(basis_to_type),
            qpoints=np.asarray(lattice_info.qpoints, dtype=float),
            num_types=self.num_types,
            output_partial=bool(self.output_partial),
            amu_2_kg=self.amu_2_kg,
        )
        self.reduced_qpoints = np.asarray(
            lattice_info.reduced_qpoints,
            dtype=float,
        )
        self._allocate_sed_accumulator(lattice_info, unique_masses)

        executor = None
        num_workers = 1
        if self.backend == "cupy":
            # CUDA is initialized only when the user selects this backend.
            self._cupy = _load_cupy()
            self._gpu_timings = {"trajectory": 0.0}
            self._gpu_event_profiler = _GPUEventProfiler(self._cupy)
            self._gpu_observed_used_bytes = 0
            self._gpu_total_bytes = 0
            self._gpu_basis_ids = tuple(
                self._cupy.asarray(atom_ids)
                for atom_ids in kernel_config.basis_ids
            )
            self._gpu_qpoints = self._cupy.asarray(
                kernel_config.qpoints
            )
            print(
                "\n****************** Using CuPy GPU backend for "
                "computing SED *****************"
            )
        elif max_cores > 1:
            # The pool lives across all blocks, and each worker retains the
            # small static config supplied by its initializer.
            # More than one worker per two Q points would create one-Q tasks,
            # so leave excess CPU cores idle for this small-Q workload.
            num_workers = min(
                max_cores,
                max(1, lattice_info.num_qpoints // 2),
            )
            executor = ProcessPoolExecutor(
                max_workers=num_workers,
                initializer=_initialize_cpu_worker,
                initargs=(kernel_config,),
            )
            print(
                f"\n****************** Using {num_workers} CPU processes "
                "for computing SED *****************"
            )

        first_block = self._get_simulation_data(
            first_trajectory_block,
            lattice_info,
        )
        remaining_blocks = (
            self._get_simulation_data(block, lattice_info)
            for block in remaining_trajectory_blocks
        )

        try:
            self._loop_over_blocks(
                lattice_info,
                kernel_config,
                executor,
                num_workers,
                chain((first_block,), remaining_blocks),
            )
        finally:
            # One shutdown replaces repeated pool startup/shutdown per block.
            if executor is not None:
                executor.shutdown()
            if self.backend == "cupy":
                # Static GPU arrays are no longer needed after the last block.
                del self._gpu_basis_ids, self._gpu_qpoints
                self._cupy.get_default_memory_pool().free_all_blocks()
                self._cupy.get_default_pinned_memory_pool().free_all_blocks()

        self._average_blocks_and_frequencies()

        elapsed = time.time() - start_time
        if self.backend == "cupy":
            self._print_gpu_timing_profile(elapsed)
            print("\n  🚀🚀 CuPy cached GPU memory has been released.")
        print(
            "\n************ Total SED calculation time: "
            f"{elapsed:.2f} seconds. ************"
        )

    def _allocate_sed_accumulator(self, lattice_info, unique_masses):
        """Allocate one block-shaped accumulator for online averaging."""

        if self.output_partial:
            print(
                "\n******* Output partial SED for each atom type and "
                "different directions *******"
            )
            print(
                "  🚀🚀 Element -> mass mapping "
                "(mass in amu; used for partial SED filenames):"
            )
            for type_index, mass in enumerate(unique_masses):
                symbol = self.type_symbols[type_index]
                print(
                    f"            {symbol} = {mass} amu"
                )

            self._sed_sum = np.zeros(
                (
                    self.num_frames_per_block,
                    lattice_info.num_qpoints,
                    self.num_types,
                    3,
                )
            )
        else:
            self._sed_sum = np.zeros(
                (
                    self.num_frames_per_block,
                    lattice_info.num_qpoints,
                )
            )

    def _loop_over_blocks(
        self,
        lattice_info,
        kernel_config,
        executor,
        num_workers,
        trajectory_blocks,
    ):
        """Read each trajectory block once and send it to the chosen backend."""

        # The discrete FFT approximates the time integral as dt * sum(...).
        scaling_const = self.dt**2 / (
            4 * np.pi * self.t_o * self.num_unit_cells
        )
        trajectory_blocks = iter(trajectory_blocks)
        for block_index in range(self.num_blocks):
            self.block_index = block_index
            self._allocate_qdot(lattice_info.num_qpoints)

            # The source supplies one exact block. With prefetch enabled, only
            # a wait for a not-yet-ready block remains on the critical path.
            read_start = time.perf_counter()
            try:
                vels, cell_vecs = next(trajectory_blocks)
            except StopIteration as error:
                raise EOFError(
                    f"Trajectory supplied only {block_index} of "
                    f"{self.num_blocks} required blocks."
                ) from error
            if getattr(self, "backend", "numpy") == "cupy":
                self._gpu_timings["trajectory"] += (
                    time.perf_counter() - read_start
                )
            if block_index == 0:
                self._print_cpu_memory_estimate(
                    vels,
                    cell_vecs,
                    kernel_config,
                    executor,
                    num_workers,
                )
                if getattr(self, "backend", "numpy") == "cupy":
                    self._check_gpu_memory(
                        vels,
                        cell_vecs,
                        kernel_config,
                    )

            print(
                "\n**************** Now calculate on averaging blocks "
                f"{block_index + 1}/{self.num_blocks} ... ****************\n"
            )
            self._loop_over_qpoints(
                vels,
                cell_vecs,
                kernel_config,
                executor,
                num_workers,
            )

            # qdot is no longer needed after this block, so scale it in place
            # and accumulate without storing every block in memory.
            self.qdot *= scaling_const
            self._sed_sum += self.qdot
            del self.qdot

    def _print_cpu_memory_estimate(
        self,
        vels,
        cell_vecs,
        kernel_config,
        executor,
        num_workers,
    ):
        """Estimate peak CPU array memory after the first block is known."""

        num_qpoints = kernel_config.qpoints.shape[0]
        batches = _qpoint_batches(num_qpoints, num_workers)
        max_q_batch = max(stop - start for start, stop in batches)
        worker_peak = _estimate_kernel_working_bytes(
            vels,
            kernel_config,
            max_q_batch,
        )

        persistent_arrays = (
            self._sed_sum.nbytes
            + self.qdot.nbytes
            + vels.nbytes
            + cell_vecs.nbytes
        )
        prefetched_block = (
            2 * vels.nbytes + cell_vecs.nbytes
            if getattr(self, "trajectory_prefetch", False)
            and getattr(self, "num_blocks", 1) > 1
            else 0
        )
        shared_velocities = vels.nbytes if executor is not None else 0
        if self.backend == "cupy":
            active_workers = 0
        else:
            active_workers = num_workers if executor is not None else 1
        estimated_bytes = (
            persistent_arrays
            + prefetched_block
            + shared_velocities
            + active_workers * worker_peak
        )
        print(
            "\n  🚀🚀 Estimated peak CPU array memory: "
            f"{estimated_bytes / 1e9:.2f} GB "
            "(trajectory block reading + SED arrays; "
            "excluding Python overhead)."
        )

    def _allocate_qdot(self, num_qpoints):
        """Allocate one block accumulator with total or xyz-resolved shape."""

        if self.output_partial:
            self.qdot = np.zeros(
                (
                    self.num_frames_per_block,
                    num_qpoints,
                    self.num_types,
                    3,
                )
            )
        else:
            self.qdot = np.zeros(
                (self.num_frames_per_block, num_qpoints)
            )

    def _loop_over_qpoints(
        self,
        vels,
        cell_vecs,
        kernel_config,
        executor,
        num_workers,
    ):
        """Dispatch a velocity block without changing its numerical meaning."""

        if self.backend == "cupy":
            self._compute_qpoints_gpu(
                vels,
                cell_vecs,
                kernel_config,
            )
        elif executor is not None:
            self._compute_qpoints_shared(
                executor,
                vels,
                cell_vecs,
                num_workers,
            )
        else:
            self._compute_qpoints_serial(
                vels,
                cell_vecs,
                kernel_config,
            )

    def _compute_qpoints_serial(self, vels, cell_vecs, kernel_config):
        """Run the batched NumPy kernel without multiprocessing."""

        num_qpoints = kernel_config.qpoints.shape[0]
        for q_start, q_stop in _qpoint_batches(num_qpoints, 1):
            self._print_qpoints(q_start, q_stop, num_qpoints)
            self.qdot[:, q_start:q_stop, ...] = _calculate_q_batch(
                vels=vels,
                cell_vecs=cell_vecs,
                qpoints=kernel_config.qpoints[q_start:q_stop],
                basis_ids_by_basis=kernel_config.basis_ids,
                config=kernel_config,
            )

    def _compute_qpoints_shared(
        self,
        executor,
        vels,
        cell_vecs,
        num_workers,
    ):
        """Share one velocity copy and let workers process small Q ranges."""

        num_qpoints = self.qdot.shape[1]
        batches = _qpoint_batches(num_qpoints, num_workers)

        # Jobs carry only a shared-memory description and Q indices, never the
        # full velocity array.
        with _share_numpy_array(vels) as shared:
            shared_name, shape, dtype_string = shared
            futures = []
            for q_start, q_stop in batches:
                self._print_qpoints(q_start, q_stop, num_qpoints)
                futures.append(
                    executor.submit(
                        _compute_shared_q_batch,
                        shared_name,
                        shape,
                        dtype_string,
                        cell_vecs,
                        q_start,
                        q_stop,
                    )
                )

            for future in as_completed(futures):
                q_start, q_stop, result = future.result()
                self.qdot[:, q_start:q_stop, ...] = result

    def _compute_qpoints_gpu(self, vels, cell_vecs, kernel_config):
        """Upload one block once and return its complete result only once."""

        cp = self._cupy
        profiler = getattr(self, "_gpu_event_profiler", None)

        upload_start = profiler.start() if profiler is not None else None
        gpu_vels = cp.asarray(vels)
        gpu_cell_vecs = cp.asarray(cell_vecs)
        if profiler is not None:
            profiler.stop("upload", upload_start)

        num_qpoints = kernel_config.qpoints.shape[0]
        batch_size = self._gpu_q_batch_size

        # Use one large batch whenever memory permits. If batching is needed,
        # retain every batch on the GPU and make only one final CPU transfer.
        if batch_size == num_qpoints:
            self._print_qpoints(0, num_qpoints, num_qpoints)
            compute_start = profiler.start() if profiler is not None else None
            gpu_result = _calculate_q_batch(
                vels=gpu_vels,
                cell_vecs=gpu_cell_vecs,
                qpoints=self._gpu_qpoints,
                basis_ids_by_basis=self._gpu_basis_ids,
                config=kernel_config,
                xp=cp,
                fft_function=cp.fft.fft,
                gpu_profiler=profiler,
            )
            if profiler is not None:
                profiler.stop("compute", compute_start)
        else:
            gpu_result = cp.empty(self.qdot.shape, dtype=cp.float64)
            for q_start in range(0, num_qpoints, batch_size):
                q_stop = min(q_start + batch_size, num_qpoints)
                self._print_qpoints(q_start, q_stop, num_qpoints)
                compute_start = (
                    profiler.start() if profiler is not None else None
                )
                batch_result = _calculate_q_batch(
                    vels=gpu_vels,
                    cell_vecs=gpu_cell_vecs,
                    qpoints=self._gpu_qpoints[q_start:q_stop],
                    basis_ids_by_basis=self._gpu_basis_ids,
                    config=kernel_config,
                    xp=cp,
                    fft_function=cp.fft.fft,
                    gpu_profiler=profiler,
                )
                gpu_result[:, q_start:q_stop, ...] = batch_result
                if profiler is not None:
                    profiler.stop("compute", compute_start)
                del batch_result

        download_start = profiler.start() if profiler is not None else None
        self.qdot[...] = cp.asnumpy(gpu_result)
        if profiler is not None:
            profiler.stop("download", download_start)

        # Release block references so CuPy can reuse these allocations.
        del gpu_vels, gpu_cell_vecs, gpu_result
        if profiler is not None:
            free_bytes, total_bytes = cp.cuda.runtime.memGetInfo()
            self._gpu_observed_used_bytes = max(
                self._gpu_observed_used_bytes,
                total_bytes - free_bytes,
            )
            self._gpu_total_bytes = total_bytes

    def _estimate_gpu_memory(
        self,
        vels,
        cell_vecs,
        kernel_config,
        q_batch_size,
    ):
        """Estimate GPU memory for one selectable Q-point batch size."""

        num_qpoints = kernel_config.qpoints.shape[0]
        working_bytes = _estimate_kernel_working_bytes(
            vels,
            kernel_config,
            q_batch_size,
        )

        # Split calculations retain the complete result on the GPU until all
        # batches finish. A single batch already is that complete result.
        retained_result = (
            self.qdot.nbytes if q_batch_size < num_qpoints else 0
        )
        return int(
            1.2
            * (
                vels.nbytes
                + cell_vecs.nbytes
                + working_bytes
                + retained_result
            )
        )

    def _check_gpu_memory(self, vels, cell_vecs, kernel_config):
        """Choose the largest safe Q batch and fail before an obvious OOM."""

        cp = self._cupy
        num_qpoints = kernel_config.qpoints.shape[0]
        free_bytes, total_bytes = cp.cuda.runtime.memGetInfo()

        # Prefer all Q points in one GPU operation. On smaller GPUs, decrease
        # only as far as required while preserving one CPU transfer per block.
        selected_batch = None
        for batch_size in range(num_qpoints, 0, -1):
            candidate_bytes = self._estimate_gpu_memory(
                vels,
                cell_vecs,
                kernel_config,
                batch_size,
            )
            if candidate_bytes <= free_bytes:
                selected_batch = batch_size
                break

        if selected_batch is None:
            raise MemoryError(
                "Estimated SED GPU memory exceeds currently available "
                "memory. Increase num_blocks or use backend = numpy."
            )

        self._gpu_q_batch_size = selected_batch
        used_bytes = total_bytes - free_bytes
        print(
            "\n  🚀🚀 GPU memory before SED: "
            f"{used_bytes / 1e9:.2f} GB used; "
            f"{free_bytes / 1e9:.2f} GB available "
            f"of {total_bytes / 1e9:.2f} GB."
        )

    def _print_gpu_timing_profile(self, total_elapsed):
        """Print one concise timing summary accumulated over all blocks."""

        event_times = self._gpu_event_profiler.elapsed_seconds()
        compute = event_times.get("compute", 0.0)
        projection = event_times.get("projection", 0.0)
        fft_time = event_times.get("fft", 0.0)
        other_gpu = max(0.0, compute - projection - fft_time)

        stages = (
            (
                "CPU: trajectory block wait",
                self._gpu_timings["trajectory"],
            ),
            ("CPU -> GPU upload", event_times.get("upload", 0.0)),
            ("GPU: tensordot projection", projection),
            ("GPU: FFT", fft_time),
            ("GPU: other operations", other_gpu),
            ("GPU -> CPU result", event_times.get("download", 0.0)),
        )
        profiled_total = sum(seconds for _, seconds in stages)
        other_sed_work = max(0.0, total_elapsed - profiled_total)

        print(
            "\n****************** CuPy SED timing profile "
            "******************"
        )
        for label, seconds in stages:
            percentage = (
                100.0 * seconds / profiled_total
                if profiled_total
                else 0.0
            )
            print(
                f"  {label:<30}: {seconds:8.2f} s "
                f"({percentage:5.1f}%)"
            )
        print(
            f"  {'Sum of the six rows above':<30}: "
            f"{profiled_total:8.2f} s"
        )
        print(
            f"  {'Other SED work':<30}: "
            f"{other_sed_work:8.2f} s"
        )
        if self._gpu_total_bytes:
            print(
                "\n  🚀🚀 Observed GPU memory after SED blocks: "
                f"{self._gpu_observed_used_bytes / 1e9:.2f} GB used "
                f"of {self._gpu_total_bytes / 1e9:.2f} GB."
            )

    def _print_qpoints(self, q_start, q_stop, num_qpoints):
        """Keep the original one-line output style for every Q point."""

        for q_index in range(q_start, q_stop):
            qpoint = self.reduced_qpoints[q_index]
            print(
                "\tNow calculating q-point "
                f"{q_index + 1}/{num_qpoints}:\t"
                f"q = ({qpoint[0]:.4f}, {qpoint[1]:.4f}, "
                f"{qpoint[2]:.4f})"
            )

    def _average_blocks_and_frequencies(self):
        """Average, fold frequencies, and convert the public SED unit."""

        # Reuse the accumulator allocation instead of creating another array.
        self._sed_sum /= self.num_blocks
        self.sed_avg = self._sed_sum
        del self._sed_sum

        n_half = len(self.freq_fft) // 2
        neg_part_flipped = self.sed_avg[:n_half:-1, ...]
        paired_length = neg_part_flipped.shape[0]
        positive = self.sed_avg[1 : 1 + paired_length, ...]
        # Fold both sides into an energy-preserving one-sided spectrum.
        # DC is kept once; paired nonzero frequencies contain +f and -f.
        self.sed_avg[1 : 1 + paired_length, ...] = (
            positive + neg_part_flipped
        )

        # Preserve the previous convention: exclude the Nyquist point.
        self.sed_avg = self.sed_avg[:n_half, ...]
        self.freq_fft = self.freq_fft[:n_half]

        # The SI calculation gives density per angular frequency in J*s.
        # Convert it to density per ordinary frequency in eV/THz so that
        # integrating the output over the THz axis directly returns eV.
        self.sed_avg *= JOULE_SECOND_TO_EV_PER_THZ

    def _get_simulation_data(self, block, lattice_info):
        """Reduce one raw trajectory block to the arrays required by SED."""

        try:
            vels = block.data[VELOCITIES]
            pos = block.data[POSITIONS]

            # Reference atoms identify corresponding cells across the block.
            cell_vecs = pos[:, lattice_info.cell_ref_ids, :].mean(axis=0)
            return vels, cell_vecs
        except Exception as error:
            raise EOFError(
                "******* Can't prepare a trajectory block\n"
                f"Error: {error}; check the trajectory and 'basis.in'. *******"
            ) from error

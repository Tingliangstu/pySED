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


"""Optional CuPy execution; importing this module does not import CuPy."""

from functools import partial

from pySED.sed.sed_kernel import calculate_q_batch
from pySED.backends.common import (
    estimate_kernel_working_bytes, print_cpu_memory_estimate, print_qpoints,
)


class _GPUEventProfiler:
    """Accumulate asynchronous CUDA event timings with minimal interference."""

    def __init__(self, cupy_module):
        """Store the CuPy runtime and initialize per-stage event records."""
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

    def measure(self, stage, operation, *args, **kwargs):
        """Enclose one queued GPU operation between CUDA events."""

        start_event = self.start()
        result = operation(*args, **kwargs)
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


def _load_cupy():
    """Import CuPy lazily and verify that a CUDA device is usable."""

    try:
        import cupy as cp
    except ImportError as error:
        raise ImportError(
            "device = gpu requires a CuPy package matching the installed "
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


class GPUBackend:
    """Compute SED blocks on one CUDA device using CuPy.

    Entering the context uploads basis indices and q points once. prepare()
    selects a q-batch size that fits the estimated free device memory.
    Each compute() uploads one trajectory block and copies its complete
    unnormalized spectrum back into the caller's NumPy result array.
    """

    def __init__(self, config, reduced_qpoints):
        """Store the topology and initialize empty device-buffer references."""
        self.config = config
        self.reduced_qpoints = reduced_qpoints
        self._cupy = None
        self._gpu_basis_ids = None
        self._gpu_qpoints = None
        self._gpu_observed_used_bytes = 0
        self._gpu_total_bytes = 0

    def __enter__(self):
        """Open the CUDA backend and upload basis indices and q points."""
        self._cupy = _load_cupy()
        try:
            self._gpu_event_profiler = _GPUEventProfiler(self._cupy)
            self._gpu_basis_ids = tuple(
                self._cupy.asarray(ids) for ids in self.config.basis_ids
            )
            self._gpu_qpoints = self._cupy.asarray(self.config.qpoints)
        except BaseException:
            self.close()
            raise
        print(
            "\n****************** Using CuPy GPU backend for "
            "computing SED *****************"
        )
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        """Release device buffers when leaving the computation context."""
        self.close()

    def close(self):
        """Release block-independent arrays and the cached device buffers."""
        self._gpu_basis_ids = None
        self._gpu_qpoints = None
        if self._cupy is not None:
            self._cupy.get_default_memory_pool().free_all_blocks()
            self._cupy.get_default_pinned_memory_pool().free_all_blocks()

    def prepare(self, vels, cell_vecs, result, accumulator, *, prefetch):
        """Report host-memory use and check available device memory."""
        print_cpu_memory_estimate(
            vels, cell_vecs, self.config, accumulator, result,
            backend="cupy", prefetch=prefetch,
        )
        self._check_gpu_memory(vels, cell_vecs, result)

    def compute(self, vels, cell_vecs, result):
        """Calculate one block on the GPU and fill the host result array in place.

        Upload velocities (frames, atoms, 3) and reference positions (cells, 3)
        once. Q batches reuse those arrays, and their results remain on the GPU
        until one final download. CUDA events record upload, projection, FFT,
        and download time without adding a synchronization after every operation.
        """

        kernel_config = self.config
        cp = self._cupy
        profiler = self._gpu_event_profiler

        project = partial(profiler.measure, "projection", cp.tensordot)
        transform = partial(profiler.measure, "fft", cp.fft.fft)
        upload_start = profiler.start()
        gpu_vels = cp.asarray(vels)
        gpu_cell_vecs = cp.asarray(cell_vecs)
        profiler.stop("upload", upload_start)

        num_qpoints = kernel_config.qpoints.shape[0]
        batch_size = self._gpu_q_batch_size

        # Use one large batch whenever memory permits. If batching is needed,
        # retain every batch on the GPU and make only one final CPU transfer.
        if batch_size == num_qpoints:
            print_qpoints(self.reduced_qpoints, 0, num_qpoints)
            compute_start = profiler.start()
            gpu_result = calculate_q_batch(
                vels=gpu_vels,
                cell_vecs=gpu_cell_vecs,
                qpoints=self._gpu_qpoints,
                basis_ids_by_basis=self._gpu_basis_ids,
                config=kernel_config,
                xp=cp,
                fft_function=transform,
                project_function=project,
            )
            profiler.stop("compute", compute_start)
        else:
            gpu_result = cp.empty(result.shape, dtype=cp.float64)
            for q_start in range(0, num_qpoints, batch_size):
                q_stop = min(q_start + batch_size, num_qpoints)
                print_qpoints(self.reduced_qpoints, q_start, q_stop)
                compute_start = profiler.start()
                batch_result = calculate_q_batch(
                    vels=gpu_vels,
                    cell_vecs=gpu_cell_vecs,
                    qpoints=self._gpu_qpoints[q_start:q_stop],
                    basis_ids_by_basis=self._gpu_basis_ids,
                    config=kernel_config,
                    xp=cp,
                    fft_function=transform,
                    project_function=project,
                )
                gpu_result[:, q_start:q_stop, ...] = batch_result
                profiler.stop("compute", compute_start)
                del batch_result

        download_start = profiler.start()
        result[...] = cp.asnumpy(gpu_result)
        profiler.stop("download", download_start)

        # Release block references so CuPy can reuse these allocations.
        del gpu_vels, gpu_cell_vecs, gpu_result
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
        result,
        q_batch_size,
    ):
        """Estimate device storage for a candidate q-batch size, in bytes.

        Include the uploaded block and kernel temporaries. When q points are
        split across batches, also count the complete result retained on the GPU.
        A 20 percent allowance covers additional allocation overhead.
        """

        kernel_config = self.config
        num_qpoints = kernel_config.qpoints.shape[0]
        working_bytes = estimate_kernel_working_bytes(
            vels,
            kernel_config,
            q_batch_size,
        )

        # Split calculations retain the complete result on the GPU until all
        # batches finish. A single batch already is that complete result.
        retained_result = (
            result.nbytes if q_batch_size < num_qpoints else 0
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

    def _check_gpu_memory(self, vels, cell_vecs, result):
        """Select the largest q batch whose estimated storage fits free GPU memory.

        Store the selected size for compute(). If even one q point does not fit,
        raise MemoryError before uploading the trajectory block. The check uses
        an allocation estimate, not a reservation of device memory.
        """

        kernel_config = self.config
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
                result,
                batch_size,
            )
            if candidate_bytes <= free_bytes:
                selected_batch = batch_size
                break

        if selected_batch is None:
            raise MemoryError(
                "Estimated SED GPU memory exceeds currently available "
                "memory. Increase num_blocks or use device = cpu."
            )

        self._gpu_q_batch_size = selected_batch
        used_bytes = total_bytes - free_bytes
        print(
            "\n  🚀🚀 GPU memory before SED: "
            f"{used_bytes / 1e9:.2f} GB used; "
            f"{free_bytes / 1e9:.2f} GB available "
            f"of {total_bytes / 1e9:.2f} GB."
        )

    def report(self, total_elapsed, trajectory_wait):
        """Print accumulated GPU-stage times alongside CPU trajectory-wait time.

        Projection and FFT timings are parts of the overall GPU compute timing;
        the remaining GPU time is reported separately to avoid counting them
        twice. All durations are seconds; memory figures use decimal GB.
        """

        event_times = self._gpu_event_profiler.elapsed_seconds()
        compute = event_times.get("compute", 0.0)
        projection = event_times.get("projection", 0.0)
        fft_time = event_times.get("fft", 0.0)
        other_gpu = max(0.0, compute - projection - fft_time)

        stages = (
            (
                "CPU: trajectory block wait",
                trajectory_wait,
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
        print("\n  🚀🚀 CuPy cached GPU memory has been released.")

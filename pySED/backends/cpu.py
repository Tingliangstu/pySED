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


"""CPU execution: serial Q batches or a persistent shared-memory pool."""

from concurrent.futures import ProcessPoolExecutor, as_completed, wait
from contextlib import contextmanager
from multiprocessing import shared_memory

import numpy as np

from pySED.sed_kernel import calculate_q_batch
from pySED.backends.common import (
    print_cpu_memory_estimate, print_qpoints, qpoint_batches,
)


_CPU_WORKER_CONFIG = None


def _initialize_cpu_worker(config):
    """Install static topology metadata once in each persistent CPU worker."""

    global _CPU_WORKER_CONFIG
    _CPU_WORKER_CONFIG = config


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
        result = calculate_q_batch(
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


class CPUBackend:
    """Own CPU workers and shared velocity buffers, never SED normalization."""

    def __init__(self, config, reduced_qpoints, max_cores):
        self.config = config
        self.reduced_qpoints = reduced_qpoints
        self.max_cores = max_cores
        self.num_workers = min(
            max_cores, max(1, config.qpoints.shape[0] // 2),
        )
        self._executor = None

    def __enter__(self):
        if self.max_cores > 1:
            self._executor = ProcessPoolExecutor(
                max_workers=self.num_workers,
                initializer=_initialize_cpu_worker,
                initargs=(self.config,),
            )
            print(
                f"\n****************** Using {self.num_workers} CPU processes "
                "for computing SED *****************"
            )
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        """One shutdown after all blocks, including a failed calculation."""
        if self._executor is not None:
            self._executor.shutdown()
            self._executor = None

    def prepare(self, vels, cell_vecs, result, accumulator, *, prefetch):
        print_cpu_memory_estimate(
            vels, cell_vecs, self.config, accumulator, result,
            backend="numpy", num_workers=self.num_workers,
            parallel=self._executor is not None, prefetch=prefetch,
        )

    def compute(self, vels, cell_vecs, result):
        """Fill one unnormalized block, serially or with persistent workers."""
        batches = qpoint_batches(len(self.config.qpoints), self.num_workers)
        if self._executor is None:
            for start, stop in batches:
                print_qpoints(self.reduced_qpoints, start, stop)
                result[:, start:stop, ...] = calculate_q_batch(
                    vels, cell_vecs, self.config.qpoints[start:stop],
                    self.config.basis_ids, self.config,
                )
            return

        # Only metadata and small Q ranges are sent through the process queue.
        with _share_numpy_array(vels) as shared:
            name, shape, dtype_string = shared
            futures = []
            try:
                for start, stop in batches:
                    print_qpoints(self.reduced_qpoints, start, stop)
                    futures.append(self._executor.submit(
                        _compute_shared_q_batch, name, shape, dtype_string,
                        cell_vecs, start, stop,
                    ))
                for future in as_completed(futures):
                    start, stop, values = future.result()
                    result[:, start:stop, ...] = values
            finally:
                # Keep the shared segment alive until every submitted reader
                # is done, even when one worker raises an exception.
                for future in futures:
                    future.cancel()
                wait(futures)

    def report(self, total_elapsed, trajectory_wait):
        """CPU has no additional device timing summary."""

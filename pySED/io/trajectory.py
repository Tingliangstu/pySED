# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
# =============================================================================

"""One block-oriented trajectory interface for SED computation."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from pathlib import Path

import numpy as np

from pySED.io.netcdf import NetCDFReader, is_netcdf
from pySED.io.schema import (
    CELLS,
    POSITIONS,
    VELOCITIES,
    TrajectoryBatch,
    TrajectoryInfo,
)
from pySED.io.text import (
    GPUMD_XYZ,
    LAMMPS_DUMP,
    iter_gpumd_xyz,
    iter_lammps_dump,
)

_SED_FIELDS = (POSITIONS, VELOCITIES, CELLS)
_END = object()


class NetCDFBlockSource:
    """Yield exact consecutive SED blocks from one NetCDF trajectory."""

    def __init__(self, path: Path, required_frames: int) -> None:
        self.path = Path(path)
        self.required_frames = int(required_frames)
        with NetCDFReader(self.path) as trajectory:
            trajectory.require(*_SED_FIELDS)
            self.info = trajectory.info
        if self.info.n_frames < self.required_frames:
            raise ValueError(
                f"Trajectory has {self.info.n_frames} frames, but "
                f"{self.required_frames} are requested."
            )

    def iter_blocks(self, block_size: int, num_blocks: int):
        """Read only the requested NetCDF frame slices, in order."""

        with NetCDFReader(self.path) as trajectory:
            trajectory.require(*_SED_FIELDS)
            for block_index in range(num_blocks):
                start = block_index * block_size
                stop = start + block_size
                data = {
                    field: trajectory.read(field, slice(start, stop))
                    for field in _SED_FIELDS
                }
                yield TrajectoryBatch(start=start, data=data)


class DirectTextBlockSource:
    """Parse one text trajectory sequentially and yield exact SED blocks."""

    def __init__(
        self,
        path: Path,
        source_format: str,
        required_frames: int,
        parser_batch_size: int,
        lammps_unit: str,
    ) -> None:
        self.path = Path(path)
        self.source_format = source_format
        self.required_frames = int(required_frames)
        self.parser_batch_size = int(parser_batch_size)
        self.lammps_unit = lammps_unit
        self.info: TrajectoryInfo | None = None
        self._started = False

    def _batches(self):
        if self.source_format == GPUMD_XYZ:
            return iter_gpumd_xyz(
                self.path,
                batch_size=self.parser_batch_size,
            )
        if self.source_format == LAMMPS_DUMP:
            return iter_lammps_dump(
                self.path,
                batch_size=self.parser_batch_size,
                lammps_unit=self.lammps_unit,
            )
        raise ValueError(
            f"direct text reading does not support '{self.source_format}'"
        )

    def iter_blocks(self, block_size: int, num_blocks: int):
        """Scan the text once and assemble the requested block boundaries."""

        if self._started:
            raise RuntimeError("a direct text trajectory can be consumed once")
        self._started = True
        expected_frames = block_size * num_blocks
        if expected_frames != self.required_frames:
            raise ValueError("trajectory block layout is inconsistent")

        batches = self._batches()
        block_data = None
        block_start = 0
        block_fill = 0
        consumed_frames = 0
        yielded_blocks = 0
        try:
            for batch in batches:
                missing = [
                    field for field in _SED_FIELDS if field not in batch.data
                ]
                if missing:
                    raise ValueError(
                        "trajectory is missing required field(s): "
                        + ", ".join(missing)
                    )
                if self.info is None:
                    positions = batch.data[POSITIONS]
                    self.info = TrajectoryInfo(
                        path=self.path,
                        n_frames=self.required_frames,
                        n_atoms=int(positions.shape[1]),
                        fields=frozenset(batch.data),
                        source_program=(
                            "GPUMD"
                            if self.source_format == GPUMD_XYZ
                            else "LAMMPS"
                        ),
                        source_format=self.source_format,
                    )

                batch_offset = 0
                while (
                    batch_offset < batch.n_frames
                    and yielded_blocks < num_blocks
                ):
                    if block_data is None:
                        block_data = {
                            field: np.empty(
                                (block_size, *batch.data[field].shape[1:]),
                                dtype=batch.data[field].dtype,
                            )
                            for field in _SED_FIELDS
                        }
                        block_start = consumed_frames
                        block_fill = 0

                    take = min(
                        block_size - block_fill,
                        batch.n_frames - batch_offset,
                    )
                    source_slice = slice(batch_offset, batch_offset + take)
                    target_slice = slice(block_fill, block_fill + take)
                    for field in _SED_FIELDS:
                        block_data[field][target_slice] = batch.data[field][
                            source_slice
                        ]
                    batch_offset += take
                    block_fill += take
                    consumed_frames += take

                    if block_fill == block_size:
                        yield TrajectoryBatch(
                            start=block_start,
                            data=block_data,
                        )
                        yielded_blocks += 1
                        block_data = None

                if yielded_blocks == num_blocks:
                    return
        finally:
            close = getattr(batches, "close", None)
            if close is not None:
                close()

        raise EOFError(
            f"Trajectory ended after {consumed_frames} requested frames; "
            f"{self.required_frames} are required."
        )


def trajectory_block_source(params):
    """Create the SED block source selected during trajectory preparation."""

    path = Path(params.trajectory_path)
    required_frames = (
        params.total_num_steps // params.output_data_stride
    )
    if is_netcdf(path):
        return NetCDFBlockSource(path, required_frames)
    return DirectTextBlockSource(
        path=path,
        source_format=params.source_format,
        required_frames=required_frames,
        parser_batch_size=params.netcdf_batch_size,
        lammps_unit=params.lammps_unit,
    )


def _next_or_end(iterator):
    try:
        return next(iterator)
    except StopIteration:
        return _END


@contextmanager
def prefetch_one(items, enabled: bool):
    """Yield an iterator with at most one item prepared in the background."""

    iterator = iter(items)
    if not enabled:
        try:
            yield iterator
        finally:
            close = getattr(iterator, "close", None)
            if close is not None:
                close()
        return

    executor = ThreadPoolExecutor(
        max_workers=1,
        thread_name_prefix="pySED-trajectory",
    )
    future = executor.submit(_next_or_end, iterator)

    def consume():
        nonlocal future
        while True:
            item = future.result()
            if item is _END:
                return
            future = executor.submit(_next_or_end, iterator)
            yield item

    consumer = consume()
    try:
        yield consumer
    finally:
        consumer.close()
        future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)
        close = getattr(iterator, "close", None)
        if close is not None:
            close()

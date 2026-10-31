# =============================================================================
#     Copyright 2025-2026 Ting Liang and pySED development team
#     This file is part of pySED.
#     pySED is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#     pySED is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#     GNU General Public License for more details.
#     You should have received a copy of the GNU General Public License
#     along with pySED. If not, see <http://www.gnu.org/licenses/>.
# =============================================================================

"""Detect and prepare the trajectory consumed by pySED."""

from __future__ import annotations

from pathlib import Path

from pySED.io.netcdf import (
    NetCDFReader,
    converted_file_is_current,
    is_netcdf,
)
from pySED.io.text import convert_text_trajectory, detect_text_format


def default_converted_path(source: Path) -> Path:
    """Return the .pysed.nc cache path beside the source trajectory."""
    path = Path(source)
    return path.with_name(f"{path.name}.pysed.nc")


def detect_trajectory(source: Path) -> str:
    """Identify text or NetCDF format from the file contents and metadata.

    NetCDF conversion metadata takes precedence over the producing program
    name. Text files are distinguished by their first header. Missing files
    raise FileNotFoundError; filenames alone do not select a parser.
    """

    source = Path(source).expanduser()
    if not source.is_file():
        raise FileNotFoundError(f"trajectory '{source}' does not exist")

    if not is_netcdf(source):
        return detect_text_format(source)

    with NetCDFReader(source) as trajectory:
        if trajectory.info.source_format:
            return trajectory.info.source_format
        program = (trajectory.info.source_program or "").lower()

    if "lammps" in program:
        return "lammps_netcdf"
    if "gpumd" in program:
        return "gpumd_netcdf"
    return "pysed_netcdf"


def prepare_trajectory(
    source: Path,
    converted_path: Path | None = None,
    batch_size: int = 32,
    compression_level: int = 1,
    storage_dtype: str = "float32",
    lammps_unit: str = "metal",
    text_mode: str = "cache",
    force: bool = False,
) -> tuple[Path, str]:
    """Return (readable_path, source_format) for the selected trajectory mode.

    Native NetCDF and direct text mode return the source path unchanged.
    Cache mode reuses a matching .pysed.nc file or rebuilds it from text when
    source metadata or conversion settings differ. force bypasses reuse only
    for text-to-NetCDF conversion. storage_dtype controls both text parsing
    and cache precision; this function does not choose the SED intensity scale.
    """

    source = Path(source).expanduser()
    source_format = detect_trajectory(source)
    if is_netcdf(source):
        return source, source_format
    if text_mode == "direct":
        return source, source_format
    if text_mode != "cache":
        raise ValueError("text_mode must be 'cache' or 'direct'")

    output = (
        Path(converted_path)
        if converted_path is not None
        else default_converted_path(source)
    )
    if force or not converted_file_is_current(
        output,
        source,
        compression_level=compression_level,
        storage_dtype=storage_dtype,
        lammps_unit=lammps_unit,
        source_format=source_format,
    ):
        convert_text_trajectory(
            source=source,
            output=output,
            source_format=source_format,
            batch_size=batch_size,
            compression_level=compression_level,
            storage_dtype=storage_dtype,
            lammps_unit=lammps_unit,
        )
    return output, source_format

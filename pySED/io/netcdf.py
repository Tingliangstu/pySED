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

"""Fast canonical NetCDF reading and streaming writing."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from pySED.io.schema import (
    ATOM_IDS,
    CELLS,
    POSITIONS,
    TIMES,
    TYPES,
    VELOCITIES,
    TrajectoryBatch,
    TrajectoryInfo,
)

_CDF_MAGICS = {b"CDF\x01", b"CDF\x02", b"CDF\x05"}
_HDF5_MAGIC = b"\x89HDF\r\n\x1a\n"
# Version 2 preserves the requested precision while parsing text trajectories.
_SCHEMA_VERSION = 2

_ALIASES = {
    POSITIONS: ("coordinates", "positions", "position"),
    VELOCITIES: ("velocities", "velocity"),
    TYPES: ("type", "types", "atom_types"),
    ATOM_IDS: ("atom_id", "atom_ids", "id"),
    TIMES: ("time", "times"),
}


def _netcdf4():
    """Import netCDF4 or raise an error with its installation command."""
    try:
        import netCDF4
    except ImportError as exc:
        raise ImportError(
            "NetCDF support requires: python -m pip install netCDF4"
        ) from exc
    return netCDF4


def is_netcdf(path: Path) -> bool:
    """Identify classic NetCDF and NetCDF4/HDF5 from file signatures."""

    try:
        with Path(path).open("rb") as stream:
            magic = stream.read(8)
    except OSError:
        return False
    return magic[:4] in _CDF_MAGICS or magic == _HDF5_MAGIC


def _first_variable(variables, names) -> str | None:
    """Return the first available variable alias, or None."""
    return next((name for name in names if name in variables), None)


def _unit(value: Any) -> str | None:
    """Normalize a unit attribute for lookup in conversion tables."""
    if value is None:
        return None
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    return (
        str(value)
        .strip()
        .lower()
        .replace("å", "angstrom")
        .replace("ångström", "angstrom")
        .replace(" ", "")
        .replace("per", "/")
    )


def _position_factor(unit: str | None) -> float:
    """Return the factor converting positions to Angstrom."""
    factors = {
        None: 1.0,
        "a": 1.0,
        "angstrom": 1.0,
        "angstroms": 1.0,
        "nm": 10.0,
        "bohr": 0.529177210903,
    }
    normalized = _unit(unit)
    if normalized not in factors:
        raise ValueError(f"unsupported position unit: {unit}")
    return factors[normalized]


def _velocity_factor(unit: str | None) -> float:
    """Return the factor converting velocities to m/s."""
    factors = {
        "m/s": 1.0,
        "angstrom/ps": 100.0,
        "angstrom/picosecond": 100.0,
        "angstrom/fs": 100000.0,
        "angstrom/femtosecond": 100000.0,
        "nm/ps": 1000.0,
    }
    normalized = _unit(unit)
    if normalized not in factors:
        raise ValueError(
            f"unsupported or missing velocity unit: {unit!r}"
        )
    return factors[normalized]


def _time_factor(unit: str | None) -> float:
    """Return the factor converting time values to ps."""
    factors = {
        None: 1.0,
        "ps": 1.0,
        "picosecond": 1.0,
        "fs": 0.001,
        "femtosecond": 0.001,
        "s": 1.0e12,
    }
    normalized = _unit(unit)
    if normalized not in factors:
        raise ValueError(f"unsupported time unit: {unit}")
    return factors[normalized]


def _cell_matrices(lengths: np.ndarray, angles: np.ndarray) -> np.ndarray:
    """Convert cell lengths and angles to row-wise lattice vectors.

    lengths and angles contain (a, b, c) and (alpha, beta, gamma) for each
    frame. Angles are in degrees; lengths keep their input unit. Return an
    array of shape (frames, 3, 3), with the first vector along x and the second
    in the xy plane. A vanishing sin(gamma) raises ValueError.
    """

    lengths = np.atleast_2d(np.asarray(lengths, dtype=float))
    angles = np.atleast_2d(np.asarray(angles, dtype=float))
    alpha, beta, gamma = np.deg2rad(angles).T
    a, b, c = lengths.T
    sin_gamma = np.sin(gamma)
    if np.any(np.isclose(sin_gamma, 0.0)):
        raise ValueError("invalid cell angle gamma")

    cells = np.zeros((len(lengths), 3, 3), dtype=float)
    cells[:, 0, 0] = a
    cells[:, 1, 0] = b * np.cos(gamma)
    cells[:, 1, 1] = b * sin_gamma
    cells[:, 2, 0] = c * np.cos(beta)
    cells[:, 2, 1] = (
        c * (np.cos(alpha) - np.cos(beta) * np.cos(gamma)) / sin_gamma
    )
    cells[:, 2, 2] = np.sqrt(
        np.maximum(
            c * c - cells[:, 2, 0] ** 2 - cells[:, 2, 1] ** 2,
            0.0,
        )
    )
    return cells


class NetCDFReader:
    """Read supported NetCDF trajectories through canonical field names.

    read() returns NumPy arrays with positions/cells in Angstrom, velocities
    in m/s, and times in ps. It reads only the requested frame slices; info
    contains metadata without loading full arrays. Use a with block or call
    close() explicitly to release the file. Reading does not reorder atoms.
    """

    def __init__(self, path: Path) -> None:
        """Open a trajectory, resolve field aliases, and validate its dimensions."""
        self.path = Path(path)
        self._dataset = _netcdf4().Dataset(str(self.path), "r")
        self._dataset.set_auto_mask(False)
        self._variables: dict[str, str] = {}

        for field, aliases in _ALIASES.items():
            name = _first_variable(self._dataset.variables, aliases)
            if name:
                self._variables[field] = name

        coordinate_name = self._variables.get(POSITIONS)
        if coordinate_name is None:
            self.close()
            raise ValueError("NetCDF trajectory has no coordinates variable")
        shape = self._dataset.variables[coordinate_name].shape
        if len(shape) != 3 or shape[-1] != 3:
            self.close()
            raise ValueError(
                f"coordinates must have shape (frame, atom, 3), found {shape}"
            )

        self._cell_vectors = _first_variable(
            self._dataset.variables,
            ("cell_vectors", "box_vectors", "box"),
        )
        self._cell_lengths = _first_variable(
            self._dataset.variables,
            ("cell_lengths",),
        )
        self._cell_angles = _first_variable(
            self._dataset.variables,
            ("cell_angles",),
        )
        fields = set(self._variables)
        if self._cell_vectors or (self._cell_lengths and self._cell_angles):
            fields.add(CELLS)

        source_program = getattr(self._dataset, "program", None)
        source_format = getattr(self._dataset, "pysed_source_format", None)
        self._info = TrajectoryInfo(
            path=self.path,
            n_frames=int(shape[0]),
            n_atoms=int(shape[1]),
            fields=frozenset(fields),
            source_program=str(source_program) if source_program else None,
            source_format=str(source_format) if source_format else None,
        )

    @property
    def info(self) -> TrajectoryInfo:
        """Return trajectory metadata without reading frame arrays."""
        return self._info

    def require(self, *fields: str) -> None:
        """Raise KeyError if any requested trajectory field is unavailable."""
        missing = sorted(set(fields) - set(self.info.fields))
        if missing:
            raise KeyError(
                f"trajectory is missing {missing}; "
                f"available fields are {sorted(self.info.fields)}"
            )

    def _slice(self, frames: slice | None) -> slice:
        """Validate a frame slice and reject reads from a closed trajectory."""
        if self._dataset is None:
            raise RuntimeError("trajectory is closed")
        if frames is None:
            return slice(None)
        if not isinstance(frames, slice):
            raise TypeError("frames must be a slice")
        return frames

    @staticmethod
    def _as_array(values, factor: float = 1.0):
        """Convert values to an array and apply a unit factor if needed."""
        array = np.asarray(values)
        return array if factor == 1.0 else array * factor

    def _read_variable(self, field: str, frames: slice):
        """Slice a frame-dependent variable or read a static variable in full."""
        variable = self._dataset.variables[self._variables[field]]
        frame_dims = {"frame", "frames", "time", "timestep", "timesteps"}
        if variable.ndim and variable.dimensions[0].lower() in frame_dims:
            return variable[frames]
        return variable[:]

    def read(self, field: str, frames: slice | None = None):
        """Read one field, optionally restricted to a slice of frames.

        field is a canonical name from io.schema. Positions/velocities retain
        shape (frames, atoms, 3), and cells use (frames, 3, 3), including a single
        frame selected by slice. Unit attributes determine conversion to Angstrom,
        m/s, or ps. Missing fields raise KeyError; unsupported units raise ValueError.
        """
        self.require(field)
        selection = self._slice(frames)
        if field == CELLS:
            return self._read_cells(selection)

        variable = self._dataset.variables[self._variables[field]]
        values = self._read_variable(field, selection)
        units = getattr(variable, "units", None)
        if field == POSITIONS:
            return self._as_array(values, _position_factor(units))
        if field == VELOCITIES:
            return self._as_array(values, _velocity_factor(units))
        if field == TIMES:
            return self._as_array(values, _time_factor(units))
        return self._as_array(values)

    def _read_cells(self, frames: slice):
        """Return selected cell matrices in Angstrom, with lattice vectors as rows.

        Repeat a static 3x3 cell for each requested frame, or read frame-dependent
        vectors directly. If only lengths and angles are stored, reconstruct
        the vectors with _cell_matrices instead.
        """
        if self._cell_vectors:
            variable = self._dataset.variables[self._cell_vectors]
            if variable.ndim == 2:
                count = len(range(*frames.indices(self.info.n_frames)))
                values = np.repeat(
                    np.asarray(variable[:])[None, :, :],
                    count,
                    axis=0,
                )
            else:
                values = variable[frames]
            return self._as_array(
                values,
                _position_factor(getattr(variable, "units", None)),
            )

        lengths_variable = self._dataset.variables[self._cell_lengths]
        angles_variable = self._dataset.variables[self._cell_angles]
        lengths = np.asarray(lengths_variable[frames])
        lengths *= _position_factor(getattr(lengths_variable, "units", None))
        cells = _cell_matrices(lengths, angles_variable[frames])
        return self._as_array(cells)

    def read_positions(self, frames: slice | None = None):
        """Read selected atom positions in Angstrom as (frames, atoms, 3)."""
        return self.read(POSITIONS, frames)

    def read_velocities(self, frames: slice | None = None):
        """Read selected atom velocities in m/s as (frames, atoms, 3)."""
        return self.read(VELOCITIES, frames)

    def read_cells(self, frames: slice | None = None):
        """Read selected cell matrices in Angstrom with lattice vectors as rows."""
        return self.read(CELLS, frames)

    def iter_batches(self, fields, batch_size: int):
        """Yield dictionaries containing the requested fields in frame batches.

        Each batch contains at most batch_size consecutive frames; the final one
        may be shorter. Values use the same canonical units as read(). Static
        fields such as atom IDs are not given an extra frame axis.
        """
        if batch_size < 1:
            raise ValueError("batch_size must be positive")
        for start in range(0, self.info.n_frames, batch_size):
            frames = slice(start, min(start + batch_size, self.info.n_frames))
            yield {field: self.read(field, frames) for field in fields}

    def close(self) -> None:
        """Close the NetCDF handle if it is open."""
        if getattr(self, "_dataset", None) is not None:
            self._dataset.close()
            self._dataset = None

    def __enter__(self):
        """Return this reader for context-managed access."""
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        """Close the trajectory when leaving the reader context."""
        self.close()


class NetCDFWriter:
    """Append trajectory batches to a NetCDF cache with explicit unit metadata.

    Input batches use the in-memory schema: Angstrom, m/s, and ps. On disk,
    velocities are stored in Angstrom/ps and converted back by NetCDFReader.
    The first append creates the file and fixes its fields and atom count.
    finish() marks successful completion; abort() leaves it incomplete.
    """

    def __init__(
        self,
        path: Path,
        source: Path,
        source_format: str | None = None,
        compression_level: int = 1,
        storage_dtype: str = "float32",
        lammps_unit: str = "metal",
    ) -> None:
        """Validate storage options and prepare an unopened trajectory writer."""
        if not 0 <= compression_level <= 9:
            raise ValueError("compression_level must be between 0 and 9")
        if storage_dtype not in {"float32", "float64"}:
            raise ValueError("storage_dtype must be float32 or float64")
        self.path = Path(path)
        self.source = Path(source)
        self.source_format = source_format
        self.compression_level = compression_level
        self.storage_dtype = storage_dtype
        self.lammps_unit = lammps_unit
        self._dataset = None
        self._variables: dict[str, Any] = {}
        self._next_frame = 0

    def _compression(self, chunks):
        """Return NetCDF chunking and optional compression settings."""
        options = {"chunksizes": chunks}
        if self.compression_level:
            options.update(
                zlib=True,
                shuffle=True,
                complevel=self.compression_level,
            )
        return options

    def _initialize(self, batch: TrajectoryBatch) -> None:
        """Create cache dimensions, fields, and provenance from the first batch.

        Record source path/size/mtime, schema version, precision, compression,
        and LAMMPS units so later runs can decide whether the cache is reusable.
        Position/velocity storage uses storage_dtype; cells and times use float64.
        The file remains marked incomplete until finish() succeeds.
        """
        positions = batch.data[POSITIONS]
        if positions.ndim != 3 or positions.shape[-1] != 3:
            raise ValueError("positions must have shape (frame, atom, 3)")
        n_atoms = int(positions.shape[1])
        float_code = "f4" if self.storage_dtype == "float32" else "f8"
        bytes_per_float = 4 if float_code == "f4" else 8
        target_bytes = 2 * 1024 * 1024
        frames_per_chunk = max(
            1,
            min(32, target_bytes // max(1, n_atoms * 3 * bytes_per_float)),
        )
        atom_chunk = max(
            1,
            min(
                n_atoms,
                target_bytes
                // max(1, frames_per_chunk * 3 * bytes_per_float),
            ),
        )
        vector_chunks = (frames_per_chunk, atom_chunk, 3)

        dataset = _netcdf4().Dataset(str(self.path), "w", format="NETCDF4")
        self._dataset = dataset
        dataset.program = "pySED"
        dataset.Conventions = "AMBER"
        dataset.pysed_schema_version = _SCHEMA_VERSION
        dataset.pysed_conversion_complete = 0
        dataset.pysed_compression_level = self.compression_level
        dataset.pysed_storage_dtype = self.storage_dtype
        dataset.pysed_lammps_unit = self.lammps_unit
        stat = self.source.resolve().stat()
        dataset.pysed_source_path = str(self.source.resolve())
        if self.source_format:
            dataset.pysed_source_format = self.source_format
        dataset.pysed_source_size = int(stat.st_size)
        dataset.pysed_source_mtime_ns = int(stat.st_mtime_ns)

        dataset.createDimension("frame", None)
        dataset.createDimension("atom", n_atoms)
        dataset.createDimension("spatial", 3)

        position_var = dataset.createVariable(
            "coordinates",
            float_code,
            ("frame", "atom", "spatial"),
            **self._compression(vector_chunks),
        )
        position_var.units = "angstrom"
        self._variables[POSITIONS] = position_var

        if VELOCITIES in batch.data:
            velocity_var = dataset.createVariable(
                "velocities",
                float_code,
                ("frame", "atom", "spatial"),
                **self._compression(vector_chunks),
            )
            velocity_var.units = "angstrom/picosecond"
            self._variables[VELOCITIES] = velocity_var

        if CELLS in batch.data:
            cell_var = dataset.createVariable(
                "cell_vectors",
                "f8",
                ("frame", "spatial", "spatial"),
                **self._compression((frames_per_chunk, 3, 3)),
            )
            cell_var.units = "angstrom"
            self._variables[CELLS] = cell_var

        if TYPES in batch.data:
            self._variables[TYPES] = dataset.createVariable(
                "type",
                "i4",
                ("frame", "atom"),
                **self._compression((frames_per_chunk, atom_chunk)),
            )

        if ATOM_IDS in batch.data:
            self._variables[ATOM_IDS] = dataset.createVariable(
                "atom_id",
                "i8",
                ("atom",),
            )

        if TIMES in batch.data:
            time_var = dataset.createVariable(
                "time",
                "f8",
                ("frame",),
                **self._compression((frames_per_chunk,)),
            )
            time_var.units = "picosecond"
            self._variables[TIMES] = time_var

    def append(self, batch: TrajectoryBatch) -> None:
        """Write the next consecutive batch, starting at the expected frame index.

        All batches must contain the same fields. Positions and cells are already
        in Angstrom; velocities are converted from m/s to the on-disk Angstrom/ps
        convention. Atom IDs are written once because their mapping is static.
        """
        if self._dataset is None:
            self._initialize(batch)
        if batch.start != self._next_frame:
            raise ValueError(
                f"expected batch starting at {self._next_frame}, got {batch.start}"
            )
        if set(batch.data) != set(self._variables):
            raise ValueError("trajectory fields changed between batches")

        frame_slice = slice(batch.start, batch.stop)
        for field, variable in self._variables.items():
            values = batch.data[field]
            if field == VELOCITIES:
                values = values / 100.0  # m/s -> Angstrom/ps
            if field == ATOM_IDS:
                if batch.start == 0:
                    variable[:] = values
            else:
                variable[frame_slice] = values
        self._next_frame = batch.stop

    def finish(self) -> None:
        """Mark a nonempty conversion complete and close its dataset."""
        if self._dataset is None or self._next_frame == 0:
            raise ValueError("trajectory contains no frames")
        self._dataset.pysed_conversion_complete = 1
        self._dataset.close()
        self._dataset = None

    def abort(self) -> None:
        """Close an unfinished dataset without marking conversion complete.

        The writer does not remove the file. The conversion function that owns
        the temporary path is responsible for removing it after an error.
        """
        if self._dataset is not None:
            self._dataset.close()
            self._dataset = None


def converted_file_is_current(
    output: Path,
    source: Path,
    compression_level: int = 1,
    storage_dtype: str = "float32",
    lammps_unit: str = "metal",
    source_format: str | None = None,
) -> bool:
    """Return whether a completed cache matches its source and conversion options.

    Compare the source path, size, and modification time plus schema version,
    storage dtype, compression, LAMMPS units, and optional source format.
    An unreadable or incomplete cache returns False. This check concerns
    trajectory conversion, not later SED settings such as legacy_unit.
    """

    if not Path(output).is_file():
        return False
    expected_path = str(Path(source).resolve())
    stat = Path(source).resolve().stat()
    try:
        with _netcdf4().Dataset(str(output), "r") as dataset:
            return (
                int(getattr(dataset, "pysed_conversion_complete", 0)) == 1
                and int(getattr(dataset, "pysed_schema_version", -1))
                == _SCHEMA_VERSION
                and getattr(dataset, "pysed_source_path", None)
                == expected_path
                and int(getattr(dataset, "pysed_source_size", -1))
                == stat.st_size
                and int(getattr(dataset, "pysed_source_mtime_ns", -1))
                == stat.st_mtime_ns
                and int(getattr(dataset, "pysed_compression_level", -1))
                == compression_level
                and getattr(dataset, "pysed_storage_dtype", None)
                == storage_dtype
                and getattr(dataset, "pysed_lammps_unit", None)
                == lammps_unit
                and (
                    source_format is None
                    or getattr(dataset, "pysed_source_format", None)
                    == source_format
                )
            )
    except (OSError, RuntimeError, ValueError):
        return False

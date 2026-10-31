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

"""Canonical trajectory field names and data containers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

POSITIONS = "positions"
VELOCITIES = "velocities"
CELLS = "cells"
TYPES = "types"
ATOM_IDS = "atom_ids"
TIMES = "times"


@dataclass(frozen=True)
class TrajectoryInfo:
    """Describe trajectory size and available fields without loading frames.

    n_frames and n_atoms are counts; fields contains the canonical names
    such as positions, velocities, and cells. Source metadata identifies the
    file format/program and is separate from the in-memory unit convention.
    """

    path: Path
    n_frames: int
    n_atoms: int
    fields: frozenset[str]
    source_program: str | None = None
    source_format: str | None = None


@dataclass
class TrajectoryBatch:
    """Hold consecutive trajectory frames using canonical field names.

    start is the zero-based index of the first frame. Positions and velocities
    have shape (frames, atoms, 3), cells (frames, 3, 3), and atom IDs (atoms,).
    Coordinates and cells are in Angstrom, velocities in m/s, and optional
    times in ps. Cell lattice vectors are rows. stop is the exclusive end index.
    """

    start: int
    data: dict[str, np.ndarray]

    @property
    def n_frames(self) -> int:
        """Return the number of frames in this batch."""
        return int(self.data[POSITIONS].shape[0])

    @property
    def stop(self) -> int:
        """Return the exclusive trajectory-frame index at the end of this batch."""
        return self.start + self.n_frames

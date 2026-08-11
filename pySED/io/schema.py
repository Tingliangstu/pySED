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
    """Trajectory metadata read without loading its arrays."""

    path: Path
    n_frames: int
    n_atoms: int
    fields: frozenset[str]
    source_program: str | None = None
    source_format: str | None = None


@dataclass
class TrajectoryBatch:
    """A consecutive block of frames produced by a text parser."""

    start: int
    data: dict[str, np.ndarray]

    @property
    def n_frames(self) -> int:
        return int(self.data[POSITIONS].shape[0])

    @property
    def stop(self) -> int:
        return self.start + self.n_frames

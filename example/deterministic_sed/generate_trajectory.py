#!/usr/bin/env python3
"""Generate a tiny deterministic LAMMPS trajectory for pySED validation."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


NUM_FRAMES = 256
DUMP_STRIDE = 10
TIME_STEP_FS = 1.0
SIGNAL_BIN = 16
SIGNAL_FREQUENCY_THz = SIGNAL_BIN / (
    NUM_FRAMES * DUMP_STRIDE * TIME_STEP_FS * 1.0e-3
)
VELOCITY_AMPLITUDE = 0.02  # Angstrom/ps for LAMMPS units metal


def atoms():
    """Return id, type, position, and the X-point phase sign."""

    result = []
    atom_id = 1
    for iz in range(2):
        for iy in range(2):
            for ix in range(2):
                result.append(
                    (atom_id, 1, ix + 0.25, iy + 0.25, iz + 0.25, (-1) ** ix)
                )
                atom_id += 1
    return result


def write_trajectory(path: Path) -> None:
    rows = atoms()
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        for frame in range(NUM_FRAMES):
            phase = 2.0 * math.pi * SIGNAL_BIN * frame / NUM_FRAMES
            # Reset the deterministic envelope in every 64-frame averaging
            # block. This broadens the line slightly while retaining its known
            # FFT-bin center in each block.
            envelope = math.exp(-(frame % 64) / 40.0)
            stream.write("ITEM: TIMESTEP\n")
            stream.write(f"{frame * DUMP_STRIDE}\n")
            stream.write("ITEM: NUMBER OF ATOMS\n8\n")
            stream.write("ITEM: BOX BOUNDS pp pp pp\n")
            stream.write("0 2\n0 2\n0 2\n")
            stream.write("ITEM: ATOMS id type x y z vx vy vz\n")
            for atom_id, atom_type, x, y, z, q_sign in rows:
                vx = VELOCITY_AMPLITUDE * q_sign * envelope * math.sin(phase)
                stream.write(
                    f"{atom_id} {atom_type} {x:.8f} {y:.8f} {z:.8f} "
                    f"{vx:.12e} 0.0 0.0\n"
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "output",
        nargs="?",
        type=Path,
        default=Path("trajectory.dump"),
        help="output path (default: trajectory.dump)",
    )
    args = parser.parse_args()
    write_trajectory(args.output)
    print(
        f"Wrote {NUM_FRAMES} frames to {args.output}; "
        f"known X-point peak = {SIGNAL_FREQUENCY_THz:.6f} THz"
    )


if __name__ == "__main__":
    main()

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


"""Compute SED from trajectory blocks using a CPU or GPU backend.

Read in this order: compute_sed -> _loop_over_blocks -> sed_kernel.
Backends only fill a block spectrum; normalization remains in the physics layer.

Trajectory readers supply positions in Angstrom and velocities in m/s.
This module prepares the crystal mapping, accumulates normalized blocks,
and exposes the final spectrum for FileIO to write.
"""

import time
from itertools import chain

import numpy as np
from scipy.fftpack import fftfreq

from pySED.backends import create_backend
from pySED.io.schema import POSITIONS, VELOCITIES
from pySED.sed.sed_kernel import (
    SEDKernelConfig, block_normalization, finalize_spectrum,
    SED_UNIT, LEGACY_SED_UNIT, LEGACY_UNIT_NOTICE,
)
from pySED.structure.atoms import atomic_masses


def _element_symbol_from_mass(mass, tolerance=5.0e-2):
    """Return the element symbol with the closest standard mass within tolerance.

    Both mass and tolerance are in amu. Return ``'unknown'`` when no match
    exists; the caller then uses a type label for partial output.
    """

    candidates = [
        (abs(mass - standard_mass), symbol)
        for symbol, standard_mass in atomic_masses.items()
        if standard_mass is not None
        and abs(mass - standard_mass) <= tolerance
    ]
    if not candidates:
        return "unknown"
    return min(candidates)[1]


class spectral_energy_density:
    """Calculate and hold the block-averaged SED for a set of q points.

    Construct with validated input parameters and trajectory metadata, then
    call ``compute_sed`` with the basis mapping and trajectory blocks.
    Results are stored in ``sed_avg`` and ``freq_fft``. Total SED has shape
    (frequency, q); partial SED adds (element/type, xyz) axes. ``sed_unit``
    identifies eV/THz or the legacy 2.x scale; frequencies are always in THz.
    """

    def __init__(self, params, trajectory_info):
        """Check atom counts and prepare the sampling interval and FFT grid.

        ``params`` supplies the time step in fs, saved-frame stride, and
        block count. ``trajectory_info`` supplies the actual atom count.
        ``dt`` and the duration ``t_o`` are stored in seconds; no trajectory
        arrays are read here. The initial FFT grid includes both signs.
        """
        self.legacy_unit = bool(getattr(params, "legacy_unit", False))
        self.sed_unit = LEGACY_SED_UNIT if self.legacy_unit else SED_UNIT

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
            f"SED output scale: {self.sed_unit} ******"
        )
        if self.legacy_unit:
            print(f"  Note: {LEGACY_UNIT_NOTICE}")

        # scipy.fftpack.fftfreq returns Hz; expose the frequency axis in THz.
        self.freq_fft = (
            fftfreq(self.num_frames_per_block, self.dt) / 1e12
        )

    def compute_sed(
        self, params, lattice_info, first_trajectory_block,
        remaining_trajectory_blocks,
    ):
        """Consume trajectory blocks and store the final spectrum on this object.

        ``lattice_info`` supplies atom-to-basis indices and Cartesian q points.
        ``first_trajectory_block`` is the block already read while preparing
        the geometry; ``remaining_trajectory_blocks`` is consumed in order.
        One backend context spans all blocks so workers or GPU buffers can
        be reused. Output files are written separately by the pipeline.
        """
        start_time = time.time()
        self.trajectory_prefetch = bool(params.trajectory_prefetch)
        kernel_config = self._prepare_kernel(params, lattice_info)

        # The trajectory source owns disk I/O. This generator only reduces
        # each supplied block to velocity and its averaged phase references.
        trajectory_blocks = (
            self._get_simulation_data(block, lattice_info)
            for block in chain((first_trajectory_block,), remaining_trajectory_blocks)
        )
        with create_backend(
            params.backend, kernel_config, self.reduced_qpoints, params.max_cores,
        ) as backend:
            self._loop_over_blocks(backend, trajectory_blocks)

        self._average_blocks_and_frequencies()
        elapsed = time.time() - start_time
        backend.report(elapsed, self._trajectory_wait)
        print(
            "\n************ Total SED calculation time: "
            f"{elapsed:.2f} seconds. ************"
        )

    def _prepare_kernel(self, params, lattice_info):
        """Build the static SEDKernelConfig and allocate the block accumulator.

        The one-based basis labels from basis.in become zero-based atom-index
        arrays, with one atom per repeated cell. Masses follow the first
        primitive-cell block, so the generated basis order must agree with
        the trajectory order. Element/type channels follow sorted unique
        masses, and Cartesian q points include the reciprocal 2*pi factor.
        """
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

        kernel_config = SEDKernelConfig(
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

        return kernel_config

    def _allocate_sed_accumulator(self, lattice_info, unique_masses):
        """Allocate the running sum on the full, two-sided FFT grid.

        Its shape is (frames_per_block, q_points), or
        (frames_per_block, q_points, types, 3) for partial output. Only one
        running sum is retained, rather than a spectrum for every block.
        """

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

    def _loop_over_blocks(self, backend, trajectory_blocks):
        """Compute each full block, apply its normalization, and add to the sum.

        ``trajectory_blocks`` yields (velocities, cell_vectors) pairs.
        The backend fills ``qdot`` in place with unnormalized FFT power;
        this method applies the physical or legacy factor before accumulation.
        Reading fewer than ``num_blocks`` blocks raises EOFError. Division
        by the block count and frequency folding happen after this loop.
        """
        scaling_const = block_normalization(
            self.dt, self.t_o, self.num_unit_cells,
            legacy_unit=getattr(self, "legacy_unit", False),
        )
        self._trajectory_wait = 0.0
        trajectory_blocks = iter(trajectory_blocks)
        for block_index in range(self.num_blocks):
            qdot = np.zeros_like(self._sed_sum)
            read_start = time.perf_counter()
            try:
                vels, cell_vecs = next(trajectory_blocks)
            except StopIteration as error:
                raise EOFError(
                    f"Trajectory supplied only {block_index} of "
                    f"{self.num_blocks} required blocks."
                ) from error
            self._trajectory_wait += time.perf_counter() - read_start

            if block_index == 0:
                backend.prepare(
                    vels, cell_vecs, qdot, self._sed_sum,
                    prefetch=self.trajectory_prefetch and self.num_blocks > 1,
                )
            print(
                "\n**************** Now calculate on averaging blocks "
                f"{block_index + 1}/{self.num_blocks} ... ****************\n"
            )
            backend.compute(vels, cell_vecs, qdot)
            qdot *= scaling_const
            self._sed_sum += qdot
            del qdot

    def _average_blocks_and_frequencies(self):
        """Finalize the accumulated spectrum and release the accumulator name.

        ``finalize_spectrum`` averages blocks, folds positive/negative pairs,
        and applies the selected output scale. The returned ``sed_avg`` is
        a view of the accumulator, so its storage remains valid after the
        private ``_sed_sum`` attribute is removed.
        """
        self.sed_avg, self.freq_fft = finalize_spectrum(
            self._sed_sum, self.freq_fft, self.num_blocks,
            legacy_unit=getattr(self, "legacy_unit", False),
        )
        del self._sed_sum

    def _get_simulation_data(self, block, lattice_info):
        """Return velocities and time-averaged cell reference positions.

        The returned arrays have shapes (frames, atoms, 3), in m/s, and
        (cells, 3), in Angstrom. ``cell_ref_ids`` selects one reference atom
        per cell; its mean position provides the spatial Fourier phase.
        Array indices must refer to the same atom order as basis.in.
        """

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

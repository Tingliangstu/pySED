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

"""Crystal geometry, basis mapping, and commensurate Q-point paths."""

# Python modules
import math
import os
from fractions import Fraction

import numpy as np

from pySED.structure import generate_data


class BZ_methods(object):
    """Build atom-to-cell metadata and the reciprocal-space path used by SED.

    Read basis.in and match the supplied primitive cell to the trajectory cell.
    After construction, qpoints contains Cartesian reciprocal vectors including
    2*pi, reduced_qpoints contains primitive reciprocal coordinates, and
    q_distances/q_labels describe the ordered path. cell_ref_ids selects one
    reference atom per unit cell for the SED spatial phase.
    """
    def __init__(
        self,
        params,
        box_in_traj,
        box_info=True,
        qpoints_info=True,
    ):

        """Read the basis mapping and match the primitive and trajectory cells."""
        self.box_info = box_info
        self.qpoints_info = qpoints_info

        # Read the unitcell and basis positions from a file
        self._get_basis_lattice(params)
        self.cell_ref_ids = np.argwhere(self.basis_index == 1)
        self.cell_ref_ids = self.cell_ref_ids.reshape(len(self.cell_ref_ids))

        # Compare lattice cell
        self._compare_cell(params, box_in_traj)

        ############## Construct Brillouin zone ################
        # calculate direct lattice vectors (for orthogonal or triclinic lattice vectors)
        self._construct_BZ_path(params)

    def _get_basis_lattice(self, params):
        """Read the four mapping columns of basis.in after its two header lines.

        Atom IDs, unit-cell indices, and basis indices are one-based in the file;
        masses are in amu. These rows must follow the trajectory atom order.
        The largest atom ID is checked against params.num_atoms before calculation.
        """
        if not os.path.exists(params.basis_lattice_file):
            raise FileNotFoundError(
                f"Basis file '{params.basis_lattice_file}' not found"
            )

        # read the lattice info from a file
        self.atom_ids, self.unitcell_index, self.basis_index, self.masses = np.loadtxt(params.basis_lattice_file,
                                                                                       skiprows=2, unpack=True)
        self.atom_ids = self.atom_ids.astype(int)
        self.unitcell_index = self.unitcell_index.astype(int)
        self.basis_index = self.basis_index.astype(int)
        self.masses = self.masses.astype(float)
        if params.num_atoms != self.atom_ids.max():
            raise ValueError("\nInconsistent input: \"num_atoms = %d\" does not match \"Max atom_id = %d\" "
                             "in basis.in file. Please carefully verify the input." %
                             (params.num_atoms, self.atom_ids.max()))

    def _reconstruct_primitive_cell(self, supercell, supercell_dim):
        """
        Reconstruct the primitive unit cell from the supercell and its dimensions.

        Parameters:
        -----------
        supercell : np.ndarray
            3x3 matrix representing the supercell.
        supercell_dim : tuple or list
            Dimensions of the supercell (number of repeats along x, y, z).

        Returns:
        --------
        primitive_cell : np.ndarray
            3x3 matrix representing the reconstructed primitive unit cell.
        """
        # Construct the scaling matrix
        scaling_matrix = np.diag(supercell_dim)

        # Inverse of the scaling matrix
        scaling_matrix_inv = np.linalg.inv(scaling_matrix)

        # Reconstruct the primitive cell
        primitive_cell = scaling_matrix_inv @ supercell

        return primitive_cell

    def _compare_cell(self, params, box_in_traj):

        ########## Output simulation cell informations ############
        """Choose the simulation cell and reconcile the primitive-cell geometry.

        Use dump_group_supercell when supplied, otherwise use the trajectory cell.
        If requested, reconstruct params.prim_unitcell from the supercell repeats,
        then apply prim_axis. These adjustments update params.prim_unitcell in place;
        subsequent reciprocal-space construction uses the resulting row vectors.
        """
        box_in_traj = np.asarray(box_in_traj, dtype=float)

        # The MD supercell can differ from the primitive cell.
        self.supercell = box_in_traj

        if getattr(params, "dump_group_supercell", None) is not None:
            if self.box_info:
                print(
                    "\n************************** Structural information "
                    "*****************************"
                )
                print("\nLocal-group SED mode is active.")
                print(
                    "   The Q-point grid is built from the local cell "
                    "specified by dump_group_supercell in the input file."
                )
                print(
                    "   Atom positions and velocities are read unchanged "
                    "from the trajectory."
                )
                print(
                    "   The full simulation cell stored in the trajectory "
                    "is not used for this local Q grid."
                )
                print(
                    "\nLocal cell specified for this dumped group "
                    "is (Angstrom):"
                )
                for row in self.supercell:
                    print(
                        "  {: .8e}    {: .8e}    {: .8e}".format(
                            row[0], row[1], row[2]
                        )
                    )
            return

        expected_prim_unitcell = params.prim_unitcell
        expected_box = generate_data.structure_maker.calculate_supercell_matrix(expected_prim_unitcell,
                                                                                params.supercell_dim)
        if params.source_format in {
            'lammps_dump',
            'lammps_netcdf',
            'gpumd_netcdf',
        }:
            expected_box = generate_data.structure_maker.calculate_restricted_cell(expected_box)
            expected_prim_unitcell = generate_data.structure_maker.calculate_restricted_cell(expected_prim_unitcell)

        if self.box_info:
            print("\n************************** Structural information *****************************")
            print("Primitive unit cell matrix is (Angstrom):")
            for row in params.prim_unitcell:
                print("  {: .8e}    {: .8e}    {: .8e}".format(row[0], row[1], row[2]))

            print("\nMD simulation cell matrix is (Angstrom):")
            for row in box_in_traj:
                print("  {: .8e}    {: .8e}    {: .8e}".format(row[0], row[1], row[2]))

        if np.allclose(box_in_traj, expected_box, atol=1e-4):
            params.prim_unitcell = expected_prim_unitcell
            if self.box_info:
                print(
                    "\nNOTE: pySED identified that structural optimization is likely performed using the NVT ensemble.")

        else:
            if self.box_info:
                print(
                    "\nWARNING: The cell matrix in the trajectory is different from calculated "
                    "by the primitive cell and supercell dimensions.")
                print("This suggests that structural optimization is likely performed using the NPT ensemble.")

            if params.rescale_prim:
                # Attempt to reconstruct the primitive cell
                params.prim_unitcell = self._reconstruct_primitive_cell(box_in_traj, params.supercell_dim)
                if self.box_info:
                    print("\nNow, pySED have reconstructed primitive unit cell matrix from MD simulation cell (Angstrom):")
                    for row in params.prim_unitcell:
                        print("  {: .8e}    {: .8e}    {: .8e}".format(row[0], row[1], row[2]))

        # see ref: https://phonopy.github.io/phonopy/setting-tags.html#primitive-axes
        if params.prim_axis is not None:
            params.prim_unitcell = params.prim_unitcell @ params.prim_axis
            if self.box_info:
                print("\nTransformation from the input unit cell to the primitive cell is performed according \'prim_axis\'.")
                print("Now primitive unit cell is transformed to (Angstrom):")
                for row in params.prim_unitcell:
                    print("  {: .8e}    {: .8e}    {: .8e}".format(row[0], row[1], row[2]))

    def _construct_BZ_path(self, params):

        # Get q_path and number of q-paths
        """Build and join the requested commensurate q-path segments.

        Input endpoints are reduced primitive reciprocal coordinates. Store
        Cartesian qpoints, reduced_qpoints, cumulative Cartesian path distances,
        and symmetry labels in matching order. Segments without allowed points
        are skipped and their labels merged into the preceding endpoint.
        """
        q_path = np.array(params.q_path, dtype=object).reshape(-1, 3)
        num_qpaths = params.num_qpaths

        try:
            BZ_Path = BZPathHelper(params.prim_unitcell, self.supercell)
        except ValueError as error:
            if getattr(params, "dump_group_supercell", None) is not None:
                raise ValueError(
                    "The local cell cannot be matched to prim_unitcell.\n"
                    "Use the local cell and primitive cell created for the "
                    "same dumped group."
                ) from error
            raise
        self._validate_dump_group_supercell(params, BZ_Path)

        qpoints_list = []
        reduced_qpoints_list = []
        q_distances = []
        q_labels = []

        label_map = {"G": r"$\Gamma$"}

        def fmt_vec(x):
            """Format a vector with four decimal places per component."""
            return "(" + ", ".join(f"{v:.4f}" for v in x) + ")"

        def shortest_distance(q1_cart, q2_cart):
            """Shortest reciprocal distance between two q-points"""
            q1_red = BZ_Path.cartesian_to_reduced(q1_cart.reshape(1, -1))[0]
            q2_red = BZ_Path.cartesian_to_reduced(q2_cart.reshape(1, -1))[0]
            dq_red = q2_red - q1_red
            dq_red -= np.round(dq_red)                # wrap to [-0.5, 0.5)
            dq_cart = BZ_Path.reduced_to_cartesian(dq_red.reshape(1, -1))[0]
            return np.linalg.norm(dq_cart)  # The Euclidean norm (length) of the difference vector

        if self.qpoints_info:
            print("\n***************************** Q-point information *****************************")

        cumulative_distance = 0.0
        last_delta = 0.0
        prev_end_label = None
        prev_end_distance = None

        for i in range(num_qpaths):

            q_start = q_path[i]
            q_end = q_path[i + 1]

            start_label = label_map.get(params.q_path_name[i], params.q_path_name[i])
            end_label = label_map.get(params.q_path_name[i + 1], params.q_path_name[i + 1])

            qpoints_cart, _ = BZ_Path.build_commensurate_path(q_start, q_end)

            if not len(qpoints_cart):
                print(f'\nWARNING: \nNo commensurate q-points found along path'
                      f' segment {i}: from {start_label}: {fmt_vec(q_start)} to {end_label}: {fmt_vec(q_end)}')

                if q_labels and prev_end_label is not None and prev_end_distance is not None:

                    merged_label = f"{prev_end_label}/{end_label}"
                    q_labels[-1] = (prev_end_distance, merged_label)

                prev_end_label = end_label
                continue

            # For distance
            distances_segment = []

            if end_label == r"$\Gamma$":
                if prev_end_label is not None and start_label != r"$\Gamma$":
                    cumulative_distance += last_delta
                if q_labels:
                    new_dist = cumulative_distance
                    old_label = q_labels[-1][1]
                    q_labels[-1] = (new_dist, old_label)

            distances_segment.append(cumulative_distance)

            for j in range(1, len(qpoints_cart)):
                step = shortest_distance(qpoints_cart[j - 1], qpoints_cart[j])
                cumulative_distance += step
                distances_segment.append(cumulative_distance)

            if len(distances_segment) >= 2:
                last_delta = distances_segment[-1] - distances_segment[-2]

            else:
                last_delta = 0.0

            # label for gamma begin
            if i == 0:
                q_labels.append((distances_segment[0], start_label))
            q_labels.append((distances_segment[-1], end_label))

            prev_end_label = end_label
            prev_end_distance = distances_segment[-1]

            # reduced coords
            reduced_qpoints = BZ_Path.cartesian_to_reduced(qpoints_cart)
            reduced_qpoints = np.where(np.isclose(reduced_qpoints, 0, atol=1e-16), 0.0, reduced_qpoints)
            qpoints_list.append(qpoints_cart)
            reduced_qpoints_list.append(reduced_qpoints)
            q_distances.extend(distances_segment)

            if self.qpoints_info:
                print(f'\nPath segment {i}: from {start_label}: {fmt_vec(q_start)} to {end_label}: {fmt_vec(q_end)}')
                print(f'Number of q-points generated: {len(reduced_qpoints)}')
                print('Reduced q-points:')
                for q in reduced_qpoints:
                    print(f'\t(reduced) q = {fmt_vec(q)}')

        self.qpoints = np.vstack(qpoints_list)
        self.reduced_qpoints = np.vstack(reduced_qpoints_list)
        self.num_qpoints = len(self.qpoints)
        self.q_distances = np.array(q_distances)
        self.q_labels = q_labels

        # End q-point information section with asterisks
        if self.qpoints_info:
            print("\n**************** The total number of q-points generated is {}  ****************".format(self.num_qpoints))

    def _validate_dump_group_supercell(self, params, bz_path):
        """Compare the local cell repetition count with the count encoded in basis.in.

        For local-group mode, abs(det(P)) must equal the largest unit-cell index.
        This is a cell-count consistency check; it does not verify every atomic
        position or atom-to-basis assignment.
        """

        if getattr(params, "dump_group_supercell", None) is None:
            return

        repetitions = abs(int(round(np.linalg.det(bz_path.P))))
        basis_cells = int(self.unitcell_index.max())
        if repetitions != basis_cells:
            raise ValueError(
                "The local cell does not match basis.in.\n"
                f"basis.in describes {basis_cells} primitive cells, but "
                f"dump_group_supercell describes {repetitions}.\n"
                "Use the local cell and basis.in created for the same "
                "dumped group."
            )


class BZPathHelper:
    """
    Helper class for handling the relationship between a primitive cell
    and a supercell, and for generating commensurate k/q-point paths in reciprocal space.

    see ref: https://gitlab.com/materials-modeling/dynasor/-/blob/master/dynasor/qpoints/lattice.py

    Features:
    - Compute an integer repetition matrix from primitive and supercell matrices
    - Calculate reciprocal lattice vectors (including the 2π factor)
    - Convert between Cartesian and reduced (fractional) reciprocal coordinates
    - Generate commensurate q-points between two reduced coordinates
    """

    def __init__(self, primitive_cell, supercell_cell):
        """Store row-wise primitive/supercell matrices and derive the repetition matrix.

        Both inputs have shape (3, 3) and use the same length unit, normally Angstrom.
        The integer matrix P satisfies P @ primitive = supercell. Raise ValueError
        when the supplied cells do not admit the required integer relation.
        """
        self._primitive_cell = np.array(primitive_cell, dtype=float)
        self._supercell_cell = np.array(supercell_cell, dtype=float)

        # Repetition matrix P such that: P @ primitive_cell = supercell_cell
        self.P = self._compute_P_matrix()

    # -------------------------
    # Properties
    # -------------------------
    @property
    def primitive(self):
        """Return the primitive cell lattice vectors as rows."""
        return self._primitive_cell

    @property
    def supercell(self):
        """Return the supercell lattice vectors as rows."""
        return self._supercell_cell

    @property
    def reciprocal_primitive(self):
        """Return reciprocal lattice of primitive cell (rows as vectors, includes 2π)."""
        return 2 * np.pi * np.linalg.inv(self.primitive.T)

    # =========================================================
    # 1. Lattice matrix computations
    # =========================================================
    def _compute_P_matrix(self):
        """Return integer P satisfying P @ primitive = supercell.

        Solve the cell relation in floating point, round the entries to integers,
        and check that rounding is within the configured numerical tolerance.
        """
        P_float = np.linalg.solve(self._primitive_cell.T, self._supercell_cell.T).T
        P_int = np.rint(P_float).astype(int)
        if not np.allclose(P_float, P_int, atol=1e-4):
            raise ValueError("Supercell is not an integer multiple of primitive cell, please check them carefully.")

        return P_int

    # -------------------------
    # Coordinate conversions
    # -------------------------
    def cartesian_to_reduced(self, qpoints_cart):
        """Convert Cartesian q-points to reduced reciprocal coordinates."""
        return np.linalg.solve(self.reciprocal_primitive.T, qpoints_cart.T).T

    def reduced_to_cartesian(self, qpoints_red):
        """Convert reduced reciprocal coordinates to Cartesian q-points."""
        return qpoints_red @ self.reciprocal_primitive

    # =========================================================
    # 3. Path generation
    # =========================================================
    def build_commensurate_path(self, start_frac, end_frac):
        """Return Cartesian q points and their fractional positions along a segment.

        start_frac/end_frac are reduced primitive reciprocal coordinates. Allowed
        positions f lie between 0 and 1; points follow start + f*(end - start).
        The returned arrays have shapes (Q, 3) and (Q,), with Cartesian reciprocal
        units including 2*pi. No allowed points produces two empty arrays.
        """
        fractions_list = self._find_allowed_fractions(start_frac, end_frac)

        if not fractions_list:
            return np.zeros((0, 3)), np.zeros((0,))

        start_float = np.asarray(start_frac, dtype=float)
        end_float = np.asarray(end_frac, dtype=float)
        q_frac_points = np.array([
            start_float + float(f) * (end_float - start_float)
            for f in fractions_list
        ])

        q_cart = self.reduced_to_cartesian(q_frac_points)
        return q_cart, np.array([float(f) for f in fractions_list])

    # =========================================================
    # 4. Internal commensurate fraction finding
    # =========================================================
    def _find_allowed_fractions(self, start_frac, end_frac):
        """Find exact fractional positions allowed by the supercell along a q segment.

        With row-vector cells satisfying P @ primitive = supercell, a reduced
        primitive reciprocal coordinate q is commensurate when q @ P.T is an
        integer vector. For q = start + f*(end - start), solve each component as
        an integer constraint on f and intersect the resulting Fraction sets.

        Return sorted fractions in [0, 1], or an empty list if the constraints
        have no solution. Coincident endpoints use the single-point convention
        and return [Fraction(0, 1)].
        """

        if np.allclose(
            np.asarray(start_frac, dtype=float),
            np.asarray(end_frac, dtype=float),
        ):
            return [Fraction(0, 1)]

        s = np.array([self._as_commensurate_fraction(x) for x in start_frac])
        e = np.array([self._as_commensurate_fraction(x) for x in end_frac])

        P_T = self.P.T
        s_sc = s @ P_T
        delta_sc = (e - s) @ P_T

        possible = None
        for a, b in zip(s_sc, delta_sc):
            results = self._solve_linear_diophantine(a, b)
            if results is None:
                continue  # no restriction along this axis
            if not results:
                return []  # no solution at all
            possible = set(results) if possible is None else possible & set(results)
        return sorted(possible) if possible else []

    def _solve_linear_diophantine(self, a, b):
        """
        Solve n = a + x*b for integers n, with a,b as Fractions and 0 <= x <= 1.

        Returns
        -------
        list(Fraction) or None
            All valid fractional x values, None if unrestricted along this axis.
        """

        if b == 0:
            return None if a.denominator == 1 else []

        # Determine min/max range for possible integer n values depending on b sign
        low, high = (a, a + b) if b > 0 else (a + b, a)

        n_min = math.floor(float(low))
        n_max = math.ceil(float(high))

        # List all integers in range and compute resulting fractions
        fracs = [Fraction(n - a, b) for n in range(n_min, n_max + 1)]

        return [f for f in fracs if 0 <= f <= 1]

    @staticmethod
    def _as_commensurate_fraction(value, max_denominator=12, atol=1e-5):
        """Convert one coordinate to the Fraction used by the q-path solver.

        Existing Fraction input is returned unchanged. For floating-point input,
        prefer a nearby small-denominator fraction within atol (for example,
        0.33333 becomes 1/3); otherwise use Fraction.limit_denominator() to form
        a rational approximation of the floating-point value.
        """
        if isinstance(value, Fraction):
            return value

        frac = Fraction(float(value))
        simple = frac.limit_denominator(max_denominator)

        if abs(float(simple) - float(value)) <= atol:
            return simple

        return frac.limit_denominator()

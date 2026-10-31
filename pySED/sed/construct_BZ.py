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
    """Generate supercell-compatible q-points on straight reciprocal paths.

    Cell vectors are rows, with ``supercell = P @ primitive``. Reduced
    q-points refer to the primitive reciprocal basis, including 2*pi in
    Cartesian coordinates. Both orthogonal and non-orthogonal cells are
    supported, including non-diagonal integer repetition matrices.

    Path intersections with the reciprocal supercell lattice are obtained
    by a parameterization based on the greatest common divisor (GCD) and
    Bézout's identity, using exact integer arithmetic. The extended Euclidean
    algorithm computes the Bézout coefficients. See
    ``_commensurate_parameters`` for the derivation.
    """

    _CELL_TOLERANCE = 1e-4

    def __init__(self, primitive_cell, supercell_cell):
        """Accept finite, nonsingular (3, 3) primitive and supercell matrices."""
        self._primitive_cell = self._validate_cell(primitive_cell, "Primitive cell")
        self._supercell_cell = self._validate_cell(supercell_cell, "Supercell")
        self.P = self._compute_P_matrix()

    @staticmethod
    def _validate_cell(cell, name):
        cell = np.array(cell, dtype=float)
        if cell.shape != (3, 3) or not np.isfinite(cell).all():
            raise ValueError(f"{name} must be a finite (3, 3) matrix.")
        if np.linalg.slogdet(cell)[0] == 0:
            raise ValueError(f"{name} must be nonsingular.")
        return cell

    @property
    def primitive(self):
        """Primitive lattice vectors as rows, in Angstrom."""
        return self._primitive_cell

    @property
    def supercell(self):
        """Supercell lattice vectors as rows, in Angstrom."""
        return self._supercell_cell

    @property
    def reciprocal_primitive(self):
        """Primitive reciprocal vectors as rows, in rad/Angstrom."""
        return 2 * np.pi * np.linalg.inv(self.primitive.T)

    def _compute_P_matrix(self):
        """Infer P in S = P C, using an absolute integer tolerance of 1e-4."""
        repetitions = np.linalg.solve(self.primitive.T, self.supercell.T).T
        nearest = np.rint(repetitions)
        if (
            not np.isfinite(nearest).all()
            or np.any(nearest < -(2**63))
            or np.any(nearest >= 2**63)
            or not np.allclose(
                repetitions, nearest, rtol=0.0, atol=self._CELL_TOLERANCE
            )
        ):
            raise ValueError(
                "Supercell must be an integer transformation of the primitive "
                "cell (absolute tolerance 1e-4)."
            )
        matrix = nearest.astype(np.int64)
        a, b, c = [[int(x) for x in row] for row in matrix]
        determinant = (
            a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0])
        )
        if determinant == 0:
            raise ValueError("The inferred integer repetition matrix is singular.")
        return matrix

    def cartesian_to_reduced(self, qpoints_cart):
        """Convert Cartesian q-points to primitive reciprocal coordinates."""
        qpoints_cart = np.asarray(qpoints_cart, dtype=float)
        return np.linalg.solve(self.reciprocal_primitive.T, qpoints_cart.T).T

    def reduced_to_cartesian(self, qpoints_red):
        """Convert primitive reciprocal coordinates to Cartesian q-points."""
        return np.asarray(qpoints_red, dtype=float) @ self.reciprocal_primitive

    def build_commensurate_path(self, start_frac, end_frac):
        """Return all commensurate points on the normalized, closed segment.

        Parameters
        ----------
        start_frac, end_frac : array-like, shape (3,)
            Endpoints in primitive reciprocal coordinates. Nearby fractions
            with denominators up to 12 are recognized within 1e-5 for decimal
            input. Explicit Fraction values are preserved exactly.

        Returns
        -------
        q_cart : ndarray, shape (N, 3)
            Cartesian q-points in rad/Angstrom, ordered from start to end.
        fractions : ndarray, shape (N,)
            Positions in [0, 1] along the normalized segment. Empty paths
            return shapes (0, 3) and (0,). A zero-length segment returns its
            point at fraction 0 only when that point is commensurate.
        """
        start = self._normalize_endpoint(start_frac)
        end = self._normalize_endpoint(end_frac)
        fractions = self._find_allowed_fractions(start, end)
        if not fractions:
            return np.zeros((0, 3)), np.zeros((0,))

        # Use the very same rational endpoints for solving and interpolation.
        # Conversion to floating point happens only after exact interpolation.
        points = np.array([
            [float(s + f * (e - s)) for s, e in zip(start, end)]
            for f in fractions
        ])
        q_cart = self.reduced_to_cartesian(points)
        turns = q_cart @ self.supercell.T / (2 * np.pi)
        if not np.isfinite(turns).all() or not np.allclose(
            turns, np.rint(turns), rtol=0.0, atol=self._CELL_TOLERANCE
        ):
            raise ValueError(
                "Generated q-points fail the supercell periodicity check "
                "(absolute tolerance 1e-4). Check the cell matrices and "
                "their numerical precision."
            )
        return q_cart, np.array([float(f) for f in fractions])

    def _find_allowed_fractions(self, start_frac, end_frac):
        """Find exact path parameters f in [0, 1] satisfying periodicity.

        With row-vector cells C and S = P @ C, a reduced q-point is
        commensurate precisely when q @ P.T is an integer vector:

            [start + f * (end - start)] @ P.T = n,  n in Z^3.

        Set a = start @ P.T and b = (end - start) @ P.T. The helper
        _commensurate_parameters solves a + f*b in Z^3 using the greatest
        common divisor (GCD) and Bézout's identity, without enumerating
        separate coordinate sets.

        Parameters
        ----------
        start_frac, end_frac : array-like, shape (3,)
            Endpoints in reduced *primitive* reciprocal coordinates.
            Decimal endpoints use the same normalization as path output.

        Returns
        -------
        list(Fraction)
            All allowed parameters, sorted in increasing order. A constant
            path returns [Fraction(0)] only if its point is commensurate.
        """
        start = self._normalize_endpoint(start_frac)
        end = self._normalize_endpoint(end_frac)
        origin = []
        direction = []
        for row in self.P:
            origin.append(sum(s * int(p) for s, p in zip(start, row)))
            direction.append(sum((e - s) * int(p) for s, e, p in zip(start, end, row)))
        return self._commensurate_parameters(origin, direction)

    @staticmethod
    def _commensurate_parameters(origin, direction):
        """Parameterize a rational line's intersection with the integer lattice.

        Let a = origin and b = direction. Clear their denominators with a
        positive common multiple D, so A = D*a and B = D*b are integers.
        For B != 0, let g = gcd(|B_i|) > 0 be the greatest common divisor
        of the components of B, and write B = g*v with a primitive integer
        direction v. The extended Euclidean algorithm computes integers w_i
        in Bézout's identity, sum(w_i*B_i) = g, hence w dot v = 1.
        Define m = w dot A.

        If n = (A + f*B)/D is integer, k = w dot n must be integer, giving

            f = (D*k - m)/g,
            n = (A - m*v)/D + k*v.

        Conversely, these equations produce every solution whenever each
        component of A - m*v is divisible by D. Otherwise the infinite
        line misses the integer lattice. The finite segment further needs

            ceil(m/D) <= k <= floor((m + g)/D).

        These bounds and the divisibility test use Python integers only.
        The output spacing is D/g; an empty integer range means no point
        on the segment. For B = 0, test A/D directly and return at most one
        representative, f = 0. No approximate endpoint-equality test is used.
        """
        denominator = math.lcm(*(x.denominator for x in (*origin, *direction)))
        offset = [int(denominator * x) for x in origin]
        step = [int(denominator * x) for x in direction]
        divisor = 0
        coefficients = [0] * len(step)
        for axis, component in enumerate(step):
            divisor, left, right = BZPathHelper._extended_gcd(divisor, component)
            coefficients = [left * x for x in coefficients]
            coefficients[axis] += right

        if divisor == 0:
            return [Fraction(0)] if all(x % denominator == 0 for x in offset) else []

        direction_int = [x // divisor for x in step]
        shift = sum(w * a for w, a in zip(coefficients, offset))
        if any((a - shift * v) % denominator for a, v in zip(offset, direction_int)):
            return []
        first = -((-shift) // denominator)
        last = (shift + divisor) // denominator
        return [Fraction(denominator * k - shift, divisor) for k in range(first, last + 1)]

    @staticmethod
    def _extended_gcd(a, b):
        """Return nonnegative g and coefficients u, v with u*a + v*b = g."""
        remainder, next_remainder = int(a), int(b)
        u, next_u, v, next_v = 1, 0, 0, 1
        while next_remainder:
            quotient = remainder // next_remainder
            remainder, next_remainder = next_remainder, remainder - quotient * next_remainder
            u, next_u = next_u, u - quotient * next_u
            v, next_v = next_v, v - quotient * next_v
        sign = 1 if remainder >= 0 else -1
        return sign * remainder, sign * u, sign * v

    @classmethod
    def _normalize_endpoint(cls, endpoint):
        values = np.asarray(endpoint, dtype=object)
        if values.shape != (3,):
            raise ValueError("Each path endpoint must contain exactly three coordinates.")
        return tuple(cls._as_commensurate_fraction(value) for value in values)

    @staticmethod
    def _as_commensurate_fraction(value, max_denominator=12, atol=1e-5):
        """Preserve exact fractions; recognize nearby simple decimal inputs."""
        try:
            frac = value if isinstance(value, Fraction) else Fraction(str(value))
            if not math.isfinite(float(frac)):
                raise ValueError
        except (ValueError, TypeError, ZeroDivisionError, OverflowError) as error:
            raise ValueError("Path coordinates must be finite real numbers or fractions.") from error
        if isinstance(value, Fraction) or (isinstance(value, str) and "/" in value):
            return frac
        simple = frac.limit_denominator(max_denominator)
        if abs(simple - frac) <= Fraction(str(atol)):
            return simple
        # Preserve pySED's Python API for values such as float(7/15), but
        # recover a fraction only when it reproduces the same float exactly.
        if isinstance(value, (float, np.floating)):
            recovered = frac.limit_denominator()
            if float(recovered) == float(value):
                return recovered
        return frac

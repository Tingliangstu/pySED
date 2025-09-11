'''
@author:
**************************  LiangTing ***************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/4/26 23:03:21 *********************
'''
# Python modules
import os
import numpy as np
import h5py
from fractions import Fraction
import math
import warnings
from pySED.structure import generate_data


class BZ_methods(object):
    def __init__(self, params, box_info=True, qpoints_info=True):

        self.box_info = box_info
        self.qpoints_info = qpoints_info

        # Read the unitcell and basis positions from a file
        self._get_basis_lattice(params)
        self.cell_ref_ids = np.argwhere(self.basis_index == 1)
        self.cell_ref_ids = self.cell_ref_ids.reshape(len(self.cell_ref_ids))

        # Compare lattice cell
        self._compare_cell(params)

        ############## Construct Brillouin zone ################
        # calculate direct lattice vectors (for orthogonal or triclinic lattice vectors)
        if params.with_eigs:  # use phonopy q-points
            pass

        else:  # create q-point list from input file (basis.in)
            # construct the BZ paths
            self._construct_BZ_path(params)

    def _get_basis_lattice(self, params):
        if not os.path.exists(params.basis_lattice_file):
            print('\nERROR: file {} not found\n'.format(params.basis_lattice_file))
            exit()

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

    def _compare_cell(self, params):

        ########## Output simulation cell informations ############
        with h5py.File(params.output_hdf5, 'r') as database:
            box_in_traj = database['box'][()]

        # Here is important, means that the supercells used for MD simulation can be different from primitive cells
        self.supercell = box_in_traj

        xhi, yhi, zhi, xy, xz, yz = generate_data.structure_maker.calculate_lattice_parameters(params.prim_unitcell,
                                                                                               params.supercell_dim)
        # Construct the expected box matrix
        expected_box = np.array([[xhi, 0.0, 0.0],
                                 [xy, yhi, 0.0],
                                 [xz, yz, zhi]])

        if self.box_info:
            print("\n************************** Structural information *****************************")
            print("Primitive unit cell matrix is (Angstrom):")
            for row in params.prim_unitcell:
                print("  {: .8e}    {: .8e}    {: .8e}".format(row[0], row[1], row[2]))

            print("\nMD simulation cell matrix is (Angstrom):")
            for row in box_in_traj:
                print("  {: .8e}    {: .8e}    {: .8e}".format(row[0], row[1], row[2]))

        if np.any(np.abs(box_in_traj - expected_box) > 1e-4):
            if self.box_info:
                print(
                    "\n⚠️Warning⚠️: The cell matrix in the trajectory is different from calculated "
                    "by the primitive cell and supercell dimensions.")
                print("This suggests that structural optimization is likely performed using the NPT ensemble.")

            if params.rescale_prim:
                # Attempt to reconstruct the primitive cell
                params.prim_unitcell = self._reconstruct_primitive_cell(box_in_traj, params.supercell_dim)
                if self.box_info:
                    print("\nNow, pySED have reconstructed primitive unit cell matrix from MD simulation cell (Angstrom):")
                    for row in params.prim_unitcell:
                        print("  {: .8e}    {: .8e}    {: .8e}".format(row[0], row[1], row[2]))
        else:
            if self.box_info:
                print(
                    "\n⚠️pySED identified that structural optimization is likely performed using the NVT ensemble.")

        # see ref: https://phonopy.github.io/phonopy/setting-tags.html#primitive-axes
        if params.prim_axis is not None:
            params.prim_unitcell = params.prim_unitcell @ params.prim_axis
            if self.box_info:
                print("\nTransformation from the input unit cell to the primitive cell is performed according \'prim_axis\'.")
                #print("If \'prim_axis\' is set, one must use NVT for MD simulation.")
                print("Now primitive unit cell is transformed to (Angstrom):")
                for row in params.prim_unitcell:
                    print("  {: .8e}    {: .8e}    {: .8e}".format(row[0], row[1], row[2]))

    def _construct_BZ_path(self, params):

        # Get q_path and number of q-paths
        q_path = np.array(params.q_path).reshape(-1, 3)
        num_qpaths = params.num_qpaths

        BZ_Path = BZPathHelper(params.prim_unitcell, self.supercell)

        qpoints_list = []
        reduced_qpoints_list = []
        q_distances = []
        q_labels = []

        label_map = {"G": "Γ"}

        def fmt_vec(x):
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
                print(f'\n⚠️Warnings⚠️: \nNo commensurate q-points found along path'
                      f' segment {i}: from {start_label}: {fmt_vec(q_start)} to {end_label}: {fmt_vec(q_end)}')

                if q_labels and prev_end_label is not None and prev_end_distance is not None:

                    #next_start_distance = None
                    # Peek directly at the global distance of the first point of the next paragraph (For test)
                    #if (i + 1) < num_qpaths:
                    #    q_start_next = q_path[i + 1]    # K
                    #    q_end_next = q_path[i + 2]      # Γ
                    #    qpoints_cart_next, _ = BZ_Path.build_commensurate_path(q_start_next, q_end_next)
                    #    if len(qpoints_cart_next) > 0:
                            # Calculate the global distance of point K
                    #        first_step_len = shortest_distance(qpoints_cart_next[0], qpoints_cart_next[1]) if len(
                    #            qpoints_cart_next) > 1 else 0.0
                    #        next_start_distance = prev_end_distance + first_step_len

                    #if next_start_distance is not None:
                    #    midpoint_distance = 0.5 * (prev_end_distance + next_start_distance)
                    #else:
                    #    midpoint_distance = prev_end_distance

                    merged_label = f"{prev_end_label}/{end_label}"
                    q_labels[-1] = (prev_end_distance, merged_label)

                    #prev_end_distance = midpoint_distance

                prev_end_label = end_label
                continue

            # For distance
            distances_segment = []

            if end_label == "Γ":
                if prev_end_label is not None and start_label != "Γ":
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
        """
        Parameters
        ----------
        primitive_cell : array-like, shape (3,3)
            Lattice vectors of the primitive cell (rows are basis vectors).
        supercell_cell : array-like, shape (3,3)
            Lattice vectors of the supercell (rows are basis vectors)
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
        """Compute integer repetition matrix R where: R @ primitive = supercell."""
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
        """
        Generate commensurate q-points along a path between two reduced coordinates.

        """
        fractions_list = self._find_allowed_fractions(start_frac, end_frac)

        if not fractions_list:
            return np.zeros((0, 3)), np.zeros((0,))

        q_frac_points = np.array([
            start_frac + float(f) * (end_frac - start_frac)
            for f in fractions_list
        ])

        q_cart = self.reduced_to_cartesian(q_frac_points)
        return q_cart, np.array([float(f) for f in fractions_list])

    # =========================================================
    # 4. Internal commensurate fraction finding
    # =========================================================
    def _find_allowed_fractions(self, start_frac, end_frac):
        """
        Find all fractional positions f in [0,1] that satisfy commensurability.
        A supercell is defined by P @ c = S for some repetition matrix P and we
        want to find fractions so that

            [start_frac + f * (end_frac - start_frac)] @ P = n

        where n is an integer multiple of the supercell size.

        Returns
        -------
        list(Fraction)
            List of allowed fractional positions f in [0,1] that satisfy commensurability.

        Parameters
        ----------
        start_frac
            start of line in reduced supercell coordinates
        end_frac
            end of line in reduced supercell coordinates
        """

        if np.allclose(start_frac, end_frac):
            return [Fraction(0, 1)]

        s = np.array([Fraction(x).limit_denominator() for x in start_frac])
        e = np.array([Fraction(x).limit_denominator() for x in end_frac])

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
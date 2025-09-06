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
            print("Primitive unit cell is (Angstrom):\n", params.prim_unitcell)
            print("\nMD simulation cell is (Angstrom):\n", box_in_traj)

        if np.any(np.abs(box_in_traj - expected_box) > 1e-4):
            if self.box_info:
                print(
                    "\nWarning: The cell in the trajectory is different from calculated by the primitive cell and supercell dimensions.")
                print("\nIf NPT ensemble is used for structural optimization, and pySED can rescale the primitive cell by setting \'rescale_prim = 1\'.")

            if params.rescale_prim:
                # Attempt to reconstruct the primitive cell
                params.prim_unitcell = self._reconstruct_primitive_cell(box_in_traj, params.supercell_dim)
                if self.box_info:
                    print("\nNow, pySED have reconstructed primitive unit cell from MD simulation cell (Angstrom):")
                    print(params.prim_unitcell)

        # https://phonopy.github.io/phonopy/setting-tags.html#primitive-axes
        if params.prim_axis is not None:

            params.prim_unitcell = params.prim_unitcell @ params.prim_axis
            if self.box_info:
                print("\nTransformation from the input unit cell to the primitive cell is performed according \'prim_axis\'.")
                #print("If \'prim_axis\' is set, one must use NVT for MD simulation.")
                print("Now primitive unit cell is transformed to (Angstrom):")
                print(params.prim_unitcell)

    def _construct_BZ_path(self, params):

        # Get q_path and number of q-paths
        q_path = np.array(params.q_path)
        num_qpaths = params.num_qpaths

        # Create Lattice instance
        lattice = Lattice(params.prim_unitcell, self.supercell)

        # Initialize lists to store qpoints
        qpoints_list = []
        reduced_qpoints_list = []
        q_segments = []                     # Will store q-points for each segment

        # Loop over each path segment
        if self.qpoints_info:
            print("\n***************************** Q-point information *****************************")

        for i in range(num_qpaths):
            q_start = q_path[i]
            q_end = q_path[i + 1]

            # Use Lattice.make_path to get commensurate qpoints along the path
            # see https://gitlab.com/materials-modeling/dynasor/-/blob/master/dynasor/qpoints/tools.py?ref_type=heads
            qpoints_cart, _ = lattice.make_path(q_start, q_end)
            if not len(qpoints_cart):
                warnings.warn(f'\nNo commensurate q-points found along path segment {i}: from {q_start} to {q_end}')
                continue

            qpoints_list.append(qpoints_cart)
            reduced_qpoints = lattice.cartesian_to_reduced(qpoints_cart)
            reduced_qpoints = np.where(np.isclose(reduced_qpoints, 0, atol=1e-16), 0.0, reduced_qpoints)
            reduced_qpoints_list.append(reduced_qpoints)
            q_segments.append(qpoints_cart)

            # Print the q-points if required
            if self.qpoints_info:
                print(f'\nPath segment {i}: from {q_start} to {q_end}')
                print(f'Number of q-points generated: {len(reduced_qpoints)}')
                print('Reduced q-points:')
                for j in range(len(reduced_qpoints)):
                    print(
                        f'\t(reduced) q = ({reduced_qpoints[j, 0]:.4f}, {reduced_qpoints[j, 1]:.4f}, {reduced_qpoints[j, 2]:.4f})')

        # Combine all qpoints
        self.qpoints = np.vstack(qpoints_list)
        self.reduced_qpoints = np.vstack(reduced_qpoints_list)
        self.num_qpoints = len(self.qpoints)                 # will be use in Phonon module

        # Calculate q_distances using cumulative distances along the path
        q_distances = [0.0]
        q_labels = dict()
        # Starting point
        qr = 0.0
        q_labels[qr] = params.q_path_name[0]        # Start label

        # Collect labels and q-distances along the entire q-point path
        idx = 0                                     # Index to track position in self.qpoints
        num_segments = len(q_segments)
        pending_label = None                        # Variable to store the label for delayed marking
        for it, q_segment in enumerate(q_segments):
            num_points_in_segment = len(q_segment)
            # Calculate the distance increment based on q-point differences
            for i in range(num_points_in_segment):
                if idx > 0:
                    delta_q = self.qpoints[idx] - self.qpoints[idx - 1]
                    qr += np.linalg.norm(delta_q)
                    q_distances.append(qr)
                    if pending_label is not None:
                        q_labels[qr] = pending_label  # Assign the label to the current distance
                        pending_label = None          # Clear the delayed label
                idx += 1
            # Add the end label of the segment
            if it < num_segments - 1:
                # For non-final segments, set the next label to the start of the next segment
                pending_label = params.q_path_name[it + 1]
            else:
                # For the final segment, mark the endpoint immediately
                q_labels[qr] = params.q_path_name[-1]

        self.q_distances = np.array(q_distances)
        self.q_labels = q_labels

        # End q-point information section with asterisks
        if self.qpoints_info:
            print("\n**************** The total number of q-points generated is {}  ****************".format(self.num_qpoints))

class Lattice(object):
    """
    The following class is borrowed from dynasor:
    https://gitlab.com/materials-modeling/dynasor/-/blob/master/dynasor/qpoints/lattice.py
    """
    def __init__(self, primitive_cell, supercell_cell):

        """Representation of a crystal supercell
        The supercell S is given by the primitive cell p and a repetition
        matrix P such that:

              dot(P, p) = S

        In this convention the cell vectors are row vectors of p and S as in
        ASE. An inverse cell is defined as:

              c_inv = inv(c).T

        and the reciprocal cell is defined as:

              c_rec = 2*pi*inv(c).T

        Notice that the inverse cell here is defined with the tranpose so that
        the lattic vectors of the inverse/reciprocal lattice are also row
        vectors. The above also implies:

              dot(P.T, S_inv) = p_inv

        The inverse supercell S_inv defines a new lattice in reciprocal space.
        Those inverse lattice points which resides inside the inverse primitive
        cell p_inv are called commensurate lattice points. These are typically
        the only points of interest in MD simulations from a crystallographic
        and lattice dynamics point of view.

        The convention taken here is that the reciprocal cell carries the 2pi
        factor onto the cartesian q-points. This is consistent with e.g.
        Kittel. The reduced coordinates are always with respect to the
        reciprocal primitive cell.

        Parameters
        ----------
        primitive
             cell metric of the primitive cell with lattice vectors as rows.
        supercell
             cell metric of the supercell with lattice vectors as rows
        """

        self._primitive_cell = np.array(primitive_cell)
        self._supercell_cell = np.array(supercell_cell)

        # Compute P matrix such that P @ primitive_cell = supercell_cell
        self.P = self._get_P_matrix()

    @property
    def primitive(self):
        """Returns the primitive cell with lattice vectors as rows"""
        return self._primitive_cell

    @property
    def supercell(self):
        """Returns the supercell with lattice vectors as rows"""
        return self._supercell_cell

    @property
    def reciprocal_primitive(self):
        """Returns inv(primitive).T so that the rows are the inverse lattice vectors"""
        return 2 * np.pi * np.linalg.inv(self.primitive.T)  # inverse lattice as rows

    @property
    def reciprocal_supercell(self):
        """Returns inv(super).T so that the rows are the inverse lattice vectors"""
        return 2 * np.pi * np.linalg.inv(self.supercell.T)  # reciprocal lattice as rows

    def _get_P_matrix(self):
        """ P c = S  ->  c.T P.T = S.T
        The P matrix must be an integer matrix
        Solve P @ primitive_cell = supercell_cell
        """
        P = np.linalg.solve(self._primitive_cell.T, self._supercell_cell.T).T
        P_rounded = np.rint(P).astype(int)
        if not np.allclose(P, P_rounded, atol=1e-4):
            raise ValueError('Supercell is not an integer multiple of primitive cell, please check them or setting \'rescale_prim = 1\'.')
        return P_rounded

    def cartesian_to_reduced(self, qpoints_cart):
        # Convert cartesian qpoints to reduced coordinates
        return np.linalg.solve(self.reciprocal_primitive.T, qpoints_cart.T).T

    def reduced_to_cartesian(self, qpoints_red):
        # Convert reduced qpoints to cartesian coordinates
        return qpoints_red @ self.reciprocal_primitive

    def make_path(self, q_start, q_end):
        """Takes qpoints in reduced coordinates and returns all points in between

        Parameters
        ----------
        start
            coordinate of starting point in reduced inverse coordinates.
            e.g. a zone mode is given as (0.5, 0, 0), (0,0,0) == (1,0,0) etc.
        stop
            stop position

        Returns
        -------
        qpoints
            coordinates of commensurate points along path in cartesian reciprocals
        dists
            fractional distance along path
        """

        fracs = self._find_on_line(q_start, q_end, self.P.T)

        if not len(fracs):
            warnings.warn('\nNo q-points along path!')
            return np.zeros((0, 3)), np.zeros((0,))
        points = np.array([q_start + float(f) * (q_end - q_start) for f in fracs])
        qpoints_cart = self.reduced_to_cartesian(points)
        dists = np.array([float(f) for f in fracs])
        return qpoints_cart, dists

    def _find_on_line(self, start, stop, P_T):
        """Find fractional distances between start and stop combatible with P

        A supercell is defined by P @ c = S for some repetition matrix P and we
        want to find fractions so that

            [start + f * (stop - start)] @ P = n

        Parameters
        ----------
        start
            start of line in reduced supercell coordinates
        stop
            end of line in reduced supercell coordinates
        P
            repetion matrix defining the supercell
        """

        if np.allclose(start, stop):
            return [Fraction(0, 1)]

        start = np.array([Fraction(s).limit_denominator() for s in start])
        stop = np.array([Fraction(s).limit_denominator() for s in stop])

        A = start @ P_T
        B = (stop - start) @ P_T

        fracs = None
        for a, b in zip(A, B):
            fs = self._solve_Diophantine(a, b)
            if fs is None:  # "inf" solutions
                continue
            elif fs == []:  # No solutions
                return []
            fracs = set(fs) if fracs is None else fracs.intersection(fs)
        return sorted(fracs)

    def _solve_Diophantine(self, a, b):
        """Solve n = a + xb for all n in Z and a,b in Q such that 0 <= x <= 1"""

        if b == 0:
            if a.denominator == 1:
                return None
            else:
                return []

        if b < 0:
            right = np.ceil(a)
            left = np.floor(a + b)
        else:
            left = np.floor(a)
            right = np.ceil(a + b)

        ns = np.arange(left, right + 1)
        fracs = [Fraction(n - a, b) for n in ns]
        fracs = [f for f in fracs if 0 <= f <= 1]

        return fracs
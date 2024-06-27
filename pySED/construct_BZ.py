'''
@author:
**************************  LiangTing ***************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/4/26 23:03:21 *********************
'''
# Python modules
import os
import numpy as np

class BZ_methods(object):
    def __init__(self, params, qpoints_info=False):
        self.params = params
        self.qpoints_info = qpoints_info
        # read the unitcell and basis positions from a file
        self._get_basis_lattice(params)
        self.cell_ref_ids = np.argwhere(self.basis_index == 1)
        self.cell_ref_ids = self.cell_ref_ids.reshape(len(self.cell_ref_ids))

        # calculate direct lattice vectors (for orthogonal or triclinic lattice vectors)
        dir_lat = params.prim_unitcell
        dir_lat[:, 0] = dir_lat[:, 0] * params.lat_params[0]
        dir_lat[:, 1] = dir_lat[:, 1] * params.lat_params[1]
        dir_lat[:, 2] = dir_lat[:, 2] * params.lat_params[2]

        # calculate reciprocal lattice vectors
        self.bz_vol = dir_lat[0, :].dot(np.cross(dir_lat[1, :], dir_lat[2, :]))
        self.recip_vecs = np.zeros((3, 3))
        self.recip_vecs[0, :] = 2 * np.pi * np.cross(dir_lat[1, :], dir_lat[2, :]) / self.bz_vol
        self.recip_vecs[1, :] = 2 * np.pi * np.cross(dir_lat[2, :], dir_lat[0, :]) / self.bz_vol
        self.recip_vecs[2, :] = 2 * np.pi * np.cross(dir_lat[0, :], dir_lat[1, :]) / self.bz_vol

        if params.with_eigs:                # use phonopy q-points
            pass
        else:                               # create q-point list from input file (basis.in)
             # construct the BZ paths
            self._construct_BZ_path(params)

        # print list of q-points to screen
        if self.qpoints_info:
            #print('\n************* There are {} q-points to be calculated *************:\n'.format(sum(params.num_qpoints)))
            for i in range(sum(params.num_qpoints)):
                print('\t    (reduced) q = ({0:.5f}, {1:.5f}, {2:.5f})'
                      .format(self.reduced_qpoints[i, 0], self.reduced_qpoints[i, 1],
                              self.reduced_qpoints[i, 2]))


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

    def _construct_BZ_path(self, params):
        # convert reduced q-points to recip. lattice-vector units
        self.q_path = np.copy(params.q_path)
        for i in range(len(params.q_path[:, 0])):
            self.q_path[i, :] = self.recip_vecs.dot(params.q_path[i, :])

        self.num_qpoints = params.num_qpoints
        self.num_qpaths = params.num_qpaths

        # Controls the user's input parameters (seems never reach here)
        if self.num_qpaths != len(self.num_qpoints):
            raise ValueError('You must specify the number of Q points on each path')

        self.qpoints = np.zeros((sum(self.num_qpoints), 3))
        self.reduced_qpoints = np.zeros((sum(self.num_qpoints), 3))

        # make q-points list by linear interpolation between symmetry points
        start = 0
        for i in range(len(self.num_qpoints)):
            end = start + self.num_qpoints[i]
            for j in range(3):
                self.qpoints[start:end, j] = np.linspace(self.q_path[i, j],
                                                         self.q_path[(i + 1), j], self.num_qpoints[i],
                                                         endpoint = False)
                self.reduced_qpoints[start:end, j] = np.linspace(params.q_path[i, j],
                                                                 params.q_path[(i + 1), j], self.num_qpoints[i],
                                                                 endpoint = False)
            start = start + self.num_qpoints[i]
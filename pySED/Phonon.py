'''
@author:
**************************  LiangTing ***************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/4/26 00:03:21 *********************
'''

import numpy as np
from scipy.fftpack import fft, fftfreq

class spectral_energy_density:
    def __init__(self, params):
        self.num_frame = params.total_num_steps // params.output_data_stride
        self.num_frame_per_split = self.num_frame // params.num_splits

        # Frequency velocities are printed (Default)
        # Its reciprocal divided by 2 is roughly the maximum frequency attainable.
        self.dt = params.time_step * params.output_data_stride / 1e15                           # from fs to second

        # total simulation time ( / 1e15 because of from fs to s)
        self.t_o = params.time_step * params.total_num_steps / params.num_splits / 1e15         # total simulation time in per splits

        # Obtaining frequencies for fourier transform
        '''
        numpy.fft.fftfreq(n, d=1.0)
        Return the Discrete Fourier Transform sample frequencies.
        f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
        f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
        '''
        self.freq_fft = fftfreq(self.num_frame_per_split, self.dt) / 1e12            # from Hz to THz

    def compute_sed(self, params, lattice_info):
        # note that velocites are in A/ps
        self.num_unit_cells = lattice_info.unitcell_index.max()
        self.num_basis = lattice_info.basis_index.max()
        self.num_loops = params.num_splits

        #do the calculation without eigenvectors
        self.sed = np.zeros((params.num_splits, self.num_frame_per_split, sum(params.num_qpoints)))
        self._loop_over_splits(params, lattice_info)
        self.sed_avg = self.sed.sum(axis=0) / self.num_loops

        max_freq = len(self.freq_fft) // 2
        self.sed_avg = self.sed_avg[:max_freq, :]
        self.freq_fft = self.freq_fft[:max_freq]

    def _loop_over_splits(self, params, lattice_info, if_scale = True):

        self.qdot = np.zeros((self.num_frame_per_split, sum(lattice_info.num_qpoints)))
        self.scale = if_scale
        if not self.scale:
            print('\n************** WARNING: You didn\'t convert any units! **************')

        for i in range(self.num_loops):
            print('\n*************** Now calculate on averaging blocks {}/{} ... ***************\n'.format(i + 1, self.num_loops))
            self.loop_index = i
            self._get_simulation_data(params, lattice_info)
            self._loop_over_qpoints(params, lattice_info)

            # unit convert (from amu to Kg)
            if if_scale:
                self.amu_2_kg = 1.66054e-27
                self.scaling_const = self.amu_2_kg / (4 * np.pi * self.t_o * self.num_unit_cells)  # J * s (unit)

            else:
                self.scaling_const = 1e-12 / (4 * np.pi * self.t_o * self.num_unit_cells)  # From s to ps

            self.sed[i, :, :] = self.qdot * self.scaling_const     # scale

    def _loop_over_qpoints(self, params, lattice_info):
        for q in range(sum(lattice_info.num_qpoints)):
            self.q_index = q
            # Output reduced_qoints for views (but use real qoints)
            print('\tNow calculate on q-point {0}/{1}:\tq = ({2:.4f}, {3:.4f}, {4:.4f})'
                  .format(q + 1, sum(lattice_info.num_qpoints), lattice_info.reduced_qpoints[q, 0],
                          lattice_info.reduced_qpoints[q, 1], lattice_info.reduced_qpoints[q, 2]))

            self.exp_fac = np.tile(lattice_info.qpoints[q, :], (self.num_unit_cells, 1))
            self.exp_fac = np.exp(1j * np.multiply(self.exp_fac, self.cell_vecs).sum(axis=1))
            self._loop_over_basis(params, lattice_info)

    def _loop_over_basis(self, params, lattice_info):

        for i in range(self.num_basis):
            basis_ids = np.argwhere(lattice_info.basis_index == (i + 1)).reshape(self.num_unit_cells)
            mass = lattice_info.masses[i]
            vx = fft(np.squeeze(self.vels[:, basis_ids, 0]) * self.exp_fac, axis = 0)
            vy = fft(np.squeeze(self.vels[:, basis_ids, 1]) * self.exp_fac, axis = 0)
            vz = fft(np.squeeze(self.vels[:, basis_ids, 2]) * self.exp_fac, axis = 0)

            # Summing over unit cells (by vx.sum(axis=1) is actually Summing over the unit cells)
            if self.scale:
                self.qdot[:, self.q_index] = (self.qdot[:, self.q_index] + (abs(vx.sum(axis=1))**2 +
                                                                           abs(vy.sum(axis=1))**2 +
                                                                           abs(vz.sum(axis=1))**2) * mass)

            else:
                self.qdot[:, self.q_index] = (self.qdot[:, self.q_index] + (abs(vx.sum(axis=1))**2 +
                                                                           abs(vy.sum(axis=1))**2 +
                                                                           abs(vz.sum(axis=1))**2) * params.time_step / 1e3
                                                                                                   * mass)    # /1e3 because of from fs to ps

    ################################################################################################
    ### read vels and pos from the hdf5 file
    def _get_simulation_data(self, params, lattice_info):
        self.vels = params.database['velocity'][self.loop_index * self.num_frame_per_split:
                                            (self.loop_index + 1) * self.num_frame_per_split, :, :]

        if self.scale:
            self.vels = self.vels * 100  # from A/ps to m/s

        self.pos = params.database['position'][self.loop_index * self.num_frame_per_split:
                                                (self.loop_index + 1) * self.num_frame_per_split, :, :]

        # time average the positions (for now, maybe can do corr. between pos and vels)
        self.cell_vecs = self.pos[:, lattice_info.cell_ref_ids, :].mean(axis = 0)

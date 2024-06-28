'''
@author:
**************************  LiangTing ***************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/4/26 00:03:21 *********************
'''

import time
import numpy as np
from scipy.fftpack import fft, fftfreq

class spectral_energy_density:
    def __init__(self, params):
        self.num_frame = params.total_num_steps // params.output_data_stride
        self.num_frame_per_split = self.num_frame // params.num_splits

        # Frequency velocities are printed (Default)
        # Its reciprocal divided by 2 is roughly the maximum frequency attainable.
        self.dt = params.time_step * params.output_data_stride / 1e15  # from fs to second

        # total simulation time ( / 1e15 because of from fs to s)
        self.t_o = params.time_step * params.total_num_steps / params.num_splits / 1e15  # total simulation time in per splits

        # unit for mass
        self.amu_2_kg = 1.66054e-27

        # for velocity conversion to generate correct SED unit (J*s)
        if params.file_format == 'lammps' and params.lammps_unit == 'metal':  ## convert velocity from A/ps to m/s (default is metal)
            self.scaling_velocity = 100
            print('\n****** Using velocity unit of A/ps in lammps (metal), the unit for SED is convert to J*s *****')

        elif params.file_format == 'lammps' and params.lammps_unit == 'real':  ## convert velocity from A/fs to m/s
            self.scaling_velocity = 100000
            print('\n****** Using velocity unit of A/fs in lammps (real), the unit for SED is convert to J*s *****')

        elif params.file_format == 'gpumd':  ## convert velocity from A/fs to m/s (for gpumd)
            self.scaling_velocity = 100000
            print('\n****** Using velocity unit of A/fs in GPUMD, the unit for SED is convert to J*s ******')

        # Obtaining frequencies for fourier transform
        '''
        numpy.fft.fftfreq(n, d=1.0)
        Return the Discrete Fourier Transform sample frequencies.
        f = [0, 1, ...,   n/2-1,     -n/2, ..., -1] / (d*n)   if n is even
        f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
        '''
        self.freq_fft = fftfreq(self.num_frame_per_split, self.dt) / 1e12  # from Hz to THz

    def compute_sed(self, params, lattice_info):

        start_time = time.time()

        # note that velocites are in A/ps
        self.num_unit_cells = lattice_info.unitcell_index.max()
        self.num_basis = lattice_info.basis_index.max()
        self.num_loops = params.num_splits

        # do the calculation without eigenvectors
        self.sed = np.zeros((params.num_splits, self.num_frame_per_split, sum(params.num_qpoints)))
        self._loop_over_splits(params, lattice_info)
        self.sed_avg = self.sed.sum(axis=0) / self.num_loops

        max_freq = len(self.freq_fft) // 2
        self.sed_avg = self.sed_avg[:max_freq, :]
        self.freq_fft = self.freq_fft[:max_freq]

        end_time = time.time()
        print(f"\n************ Time for SED computing taken: {end_time - start_time} seconds. ************")

    def _loop_over_splits(self, params, lattice_info, if_get_memory=True):

        if if_get_memory:
            import psutil
            process = psutil.Process()
            initial_memory = process.memory_info().rss
            single_thread_memory_usage = initial_memory / 1e6
            estimated_memory_usage = single_thread_memory_usage
            print('\n**************** Estimated computing memory usage: {:.2f} MB. ***************'.format(estimated_memory_usage))

        self.qdot = np.zeros((self.num_frame_per_split, sum(lattice_info.num_qpoints)))

        for i in range(self.num_loops):
            print('\n**************** Now calculate on averaging blocks {}/{} ... ****************\n'.format(i + 1,
                                                                                                        self.num_loops))
            self.loop_index = i
            self._get_simulation_data(params, lattice_info)
            self._loop_over_qpoints(lattice_info)

            self.scaling_const = 1 / (4 * np.pi * self.t_o * self.num_unit_cells)
            # self.qdot in kg*m^2, then /self.t_o, so convert to kg*m^2/s, and finally convert to J * s (J=m^2·kg·s-2)
            self.sed[i, :, :] = self.qdot * self.scaling_const  # scale

    def _loop_over_qpoints(self, lattice_info):
        for q in range(sum(lattice_info.num_qpoints)):
            self.q_index = q
            # Output reduced_qoints for views (but use real qoints)
            print('\tNow calculate on q-point {0}/{1}:\tq = ({2:.4f}, {3:.4f}, {4:.4f})'
                  .format(q + 1, sum(lattice_info.num_qpoints), lattice_info.reduced_qpoints[q, 0],
                          lattice_info.reduced_qpoints[q, 1], lattice_info.reduced_qpoints[q, 2]))

            self.exp_fac = np.tile(lattice_info.qpoints[q, :], (self.num_unit_cells, 1))
            self.exp_fac = np.exp(1j * np.multiply(self.exp_fac, self.cell_vecs).sum(axis=1))
            self._loop_over_basis(lattice_info)

    def _loop_over_basis(self, lattice_info):
        for i in range(self.num_basis):
            basis_ids = np.argwhere(lattice_info.basis_index == (i + 1)).reshape(self.num_unit_cells)
            mass = lattice_info.masses[i]
            vx = fft(np.squeeze(self.vels[:, basis_ids, 0]) * self.exp_fac, axis=0)
            vy = fft(np.squeeze(self.vels[:, basis_ids, 1]) * self.exp_fac, axis=0)
            vz = fft(np.squeeze(self.vels[:, basis_ids, 2]) * self.exp_fac, axis=0)
            # Sum over all the unit cells after the FFT, before computing the modulus squared
            vx_sum = vx.sum(axis=1)
            vy_sum = vy.sum(axis=1)
            vz_sum = vz.sum(axis=1)

            # Compute the sum of the squares of the modulus of the summed FFT components
            mod_squared = np.abs(vx_sum) ** 2 + np.abs(vy_sum) ** 2 + np.abs(vz_sum) ** 2

            # Scale the result by the mass of the basis and update the spectral energy density (qdot)
            self.qdot[:, self.q_index] += mod_squared * mass * self.amu_2_kg  # Unit conversion (from amu to Kg)

    ################################################################################################
    ### read vels and pos from the hdf5 file

    def _get_simulation_data(self, params, lattice_info):

        self.vels = params.database['velocity'][self.loop_index * self.num_frame_per_split:
                                                (self.loop_index + 1) * self.num_frame_per_split, :, :]

        self.vels = self.vels * self.scaling_velocity  # from A/ps to m/s or from A/fs to m/s

        self.pos = params.database['position'][self.loop_index * self.num_frame_per_split:
                                               (self.loop_index + 1) * self.num_frame_per_split, :, :]

        # time average the positions (for now, maybe can do corr. between pos and vels)
        self.cell_vecs = self.pos[:, lattice_info.cell_ref_ids, :].mean(axis=0)

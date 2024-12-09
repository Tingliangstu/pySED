'''
@author:
**************************  LiangTing ***************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/4/26 00:03:21 *********************
'''

import time
import numpy as np
from scipy.fftpack import fft, fftfreq
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import h5py

class spectral_energy_density:
    def __init__(self, params):
        self.num_frame = params.total_num_steps // params.output_data_stride
        self.num_frame_per_split = self.num_frame // params.num_splits

        print('\nThe number of frames used to calculate SED is {0}.'.format(self.num_frame))
        print('\nThe {0} frame trajectories were divided into {1} blocks for averaging.'.format(self.num_frame, params.num_splits))

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

        # For multi-threading
        self.use_parallel = params.use_parallel
        self.max_cores = params.max_cores
        if self.use_parallel and self.max_cores is None:
            self.max_cores = os.cpu_count()

        if self.use_parallel and self.max_cores > 1:
            print('\n****************** Using {0} cores parallelism for computing SED *****************'.format(
                self.max_cores))

        start_time = time.time()

        # note that velocites are in A/ps
        self.num_unit_cells = lattice_info.unitcell_index.max()
        self.num_basis = lattice_info.basis_index.max()
        self.num_loops = params.num_splits

        # do the calculation without eigenvectors
        self.sed = np.zeros((params.num_splits, self.num_frame_per_split, lattice_info.num_qpoints))
        self._loop_over_splits(params, lattice_info)
        self.sed_avg = self.sed.sum(axis=0) / self.num_loops

        max_freq = len(self.freq_fft) // 2
        self.sed_avg = self.sed_avg[:max_freq, :]
        self.freq_fft = self.freq_fft[:max_freq]

        end_time = time.time()
        print(f"\n************ Time for SED computing taken: {end_time - start_time:.2f} seconds. ************")

    def _loop_over_splits(self, params, lattice_info, if_get_memory=True):

        if if_get_memory:
            import psutil
            process = psutil.Process()
            initial_memory = process.memory_info().rss
            single_thread_memory_usage = initial_memory / 1e6
            estimated_memory_usage = single_thread_memory_usage * (self.max_cores if self.use_parallel else 1)
            print('\n**************** Estimated computing memory usage: {:.2f} MB. ***************'.format(estimated_memory_usage))

        for i in range(self.num_loops):
            print('\n**************** Now calculate on averaging blocks {}/{} ... ****************\n'.format(i + 1,
                                                                                                        self.num_loops))
            self.loop_index = i
            self.qdot = np.zeros((self.num_frame_per_split, lattice_info.num_qpoints))
            vels, cell_vecs = self._get_simulation_data(params, lattice_info)
            self._loop_over_qpoints(lattice_info, vels, cell_vecs)

            self.scaling_const = 1 / (4 * np.pi * self.t_o * self.num_unit_cells)
            # self.qdot in kg*m^2, then /self.t_o, so convert to kg*m^2/s, and finally convert to J * s (J=m^2·kg·s-2)
            self.sed[i, :, :] = self.qdot * self.scaling_const  # scale

    def _loop_over_qpoints(self, lattice_info, vels, cell_vecs):

        q_indices = list(range(lattice_info.num_qpoints))

        if self.use_parallel and self.max_cores > 1:
            args_list = [(q, lattice_info, vels, cell_vecs) for q in q_indices]
            with ProcessPoolExecutor(max_workers=self.max_cores) as executor:
                futures = [executor.submit(self.process_q_point, *args) for args in args_list]
                for future in as_completed(futures):
                    q_index, qdot_q = future.result()
                    self.qdot[:, q_index] = qdot_q
        else:                                             # use one core for windows system
            for q in q_indices:
                q_index, qdot_q = self.process_q_point(q, lattice_info, vels, cell_vecs)
                self.qdot[:, q_index] = qdot_q

    def process_q_point(self, q_index, lattice_info, vels, cell_vecs):

        print('\tNow calculating q-point {0}/{1}:\tq = ({2:.4f}, {3:.4f}, {4:.4f})'
              .format(q_index + 1, lattice_info.num_qpoints,
                      lattice_info.reduced_qpoints[q_index, 0],
                      lattice_info.reduced_qpoints[q_index, 1],
                      lattice_info.reduced_qpoints[q_index, 2]))

        #exp_fac = np.tile(lattice_info.qpoints[q_index, :], (self.num_unit_cells, 1))
        #exp_fac = np.exp(1j * np.multiply(exp_fac, cell_vecs).sum(axis=1))
        exp_fac = np.exp(1.0j * np.dot(cell_vecs, lattice_info.qpoints[q_index, :]))   # for triclinic cell

        qdot_q = self._loop_over_basis(vels, exp_fac, lattice_info)

        return q_index, qdot_q

    def _loop_over_basis(self, vels, exp_fac, lattice_info):

        qdot_q = np.zeros(vels.shape[0])

        for i in range(self.num_basis):
            basis_ids = np.argwhere(lattice_info.basis_index == (i + 1)).reshape(self.num_unit_cells)
            mass = lattice_info.masses[i]
            vx = fft(np.squeeze(vels[:, basis_ids, 0]) * exp_fac, axis=0)
            vy = fft(np.squeeze(vels[:, basis_ids, 1]) * exp_fac, axis=0)
            vz = fft(np.squeeze(vels[:, basis_ids, 2]) * exp_fac, axis=0)
            # Sum over all the unit cells after the FFT, before computing the modulus squared
            vx_sum = vx.sum(axis=1)
            vy_sum = vy.sum(axis=1)
            vz_sum = vz.sum(axis=1)

            # Compute the sum of the squares of the modulus of the summed FFT components
            mod_squared = np.abs(vx_sum) ** 2 + np.abs(vy_sum) ** 2 + np.abs(vz_sum) ** 2

            # Scale the result by the mass of the basis and update the spectral energy density (qdot-sum over basis)
            qdot_q += mod_squared * mass * self.amu_2_kg  # Unit conversion (from amu to Kg)
        return qdot_q

    ################################################################################################
    ### read vels and pos from the hdf5 file

    def _get_simulation_data(self, params, lattice_info):
        try:
            # Open the HDF5 file in read-only mode
            with h5py.File(params.output_hdf5, 'r') as database:
                vels = database['velocity'][self.loop_index * self.num_frame_per_split:
                                            (self.loop_index + 1) * self.num_frame_per_split, :, :]

                vels = vels * self.scaling_velocity  # From A/ps to m/s or from A/fs to m/s

                pos = database['position'][self.loop_index * self.num_frame_per_split:
                                           (self.loop_index + 1) * self.num_frame_per_split, :, :]

            # Time average the positions
            cell_vecs = pos[:, lattice_info.cell_ref_ids, :].mean(axis=0)

            return vels, cell_vecs

        except Exception as e:
            raise EOFError(f'******* Can\'t open {params.output_hdf5}\nError: {e}, please check its integrity and that of \'basis.in\' file. *******')


'''
@author:
**************************  LiangTing ***************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/4/26 01:03:21 *********************
'''

import numpy as np
import os
import h5py

def print_error(txt):
    print('\nERROR: Value for input paramater {} seems wrong or not exist, see README.\n'.format(txt))
    exit()

class get_parse_input(object):
    def __init__(self, input_file = 'input_SED.in'):

        self.input_file = input_file
        ## ************ Control parameters ************
        self.compress = False
        self.num_splits = 1

        ### ************ MD simulation parameters **************
        self.num_atoms = 0
        self.total_num_steps = 0
        self.time_step = 0
        self.output_data_stride = 0

        ### ************ Input and output files **************
        self.vels_file = 'vels.dat'
        self.pos_file = 'pos.dat'
        self.basis_lattice_file = 'basis.in'
        self.output_hdf5 = 'vel_pos_compress.hdf5'
        self.file_format = 'lammps'

        self.plot_lorentz = False
        self.plot_cutoff_freq = None
        self.plot_interval = 5            # Thz
        self.peak_max_hwhm = 1e6          # default

        ## ********** eigenvector from phonopy for further development
        self.with_eigs = None

        input_txt = open(self.input_file, 'r').readlines()

        for line in input_txt:
            # skip blank and comment lines
            if len(line.strip()) == 0:
                continue
            elif line.strip()[0] == '#':
                continue
            # control parameters
            txt = line.strip().split()

            # MD simulation control parameters
            if txt[0] == 'num_atoms':
                try:
                    self.num_atoms = int(txt[txt.index('=') + 1])    # The index of label('=') plus 1
                except:
                    print_error('num_atoms')
            elif txt[0] == 'total_num_steps':
                try:
                    self.total_num_steps = int(txt[txt.index('=') + 1])  # The index of label('=') plus 1
                except:
                    print_error('total_num_steps')
            elif txt[0] == 'time_step':
                try:
                    self.time_step = float(txt[txt.index('=') + 1])  # The index of label('=') plus 1
                except:
                    print_error('time_step')
            elif txt[0] == 'output_data_stride':
                try:
                    self.output_data_stride = int(txt[txt.index('=') + 1])  # The index of label('=') plus 1
                except:
                    print_error('output_data_stride')

            # Lattice dynamic
            elif txt[0] == 'lat_params':
                try:
                    self.lat_params = np.array(txt[(txt.index('=') + 1):
                                                   (txt.index('=') + 4)]).astype(float)
                except:
                    print_error('lat_params')
            elif txt[0] == 'prim_unitcell':
                try:
                    self.prim_unitcell = np.array(txt[(txt.index('=') + 1):
                                                  (txt.index('=') + 10)]).astype(float)
                    self.prim_unitcell = self.prim_unitcell.reshape(3, 3)

                except:
                    print_error('prim_unitcell')

            # Q-points
            elif txt[0] == 'num_qpaths':
                try:
                    self.num_qpaths = int(txt[txt.index('=')+1])
                except:
                    print_error('num_qpaths')

            elif txt[0] == 'num_qpoints':
                try:
                    self.num_qpoints = np.array(txt[(txt.index('=')+1):(txt.index('=')+self.num_qpaths+1)]).astype(int)
                except:
                    print_error('num_qpoints')

            elif txt[0] == 'q_path':
                try:
                    self.q_path = np.array(txt[(txt.index('=')+1):
                        (txt.index('=')+int((self.num_qpaths+1)*3+1))]).astype(float)
                    self.q_path = self.q_path.reshape(self.num_qpaths+1, 3)

                except:
                    print_error('q_path')

            # whether or not to build hdf5 database
            elif txt[0] == 'compress':
                try:
                    self.compress = bool(int(txt[txt.index('=') + 1]))
                except:
                    print_error('compress')

            # number of blocks to average.
            elif txt[0] == 'num_splits':
                try:
                    self.num_splits = int(txt[txt.index('=') + 1])
                except:
                    print_error('num_splits')

            # change the file name for input and output
            elif txt[0] == 'basis_lattice_file':
                try:
                    self.basis_lattice_file = str(txt[txt.index('=') + 1].strip('\''))
                except:
                    print_error('basis_lattice_file')
            elif txt[0] == 'vels_file':
                try:
                    self.vels_file = str(txt[txt.index('=') + 1].strip('\''))
                except:
                    print_error('vels_file')

            elif txt[0] == 'pos_file':
                try:
                    self.pos_file = str(txt[txt.index('=') + 1].strip('\''))
                except:
                    print_error('pos_file')

            elif txt[0] == 'output_hdf5':
                try:
                    self.output_hdf5 = str(txt[txt.index('=') + 1].strip('\''))
                except:
                    print_error('output_hdf5')

            elif txt[0] == 'out_files_name':
                try:
                    self.out_files_name = str(txt[txt.index('=')+1].strip('\''))
                except:
                    print_error('out_files_name')

            elif txt[0] == 'file_format':
                try:
                    self.file_format = str(txt[txt.index('=') + 1].strip('\''))
                except:
                    print_error('file_format')

            # plotting
            elif txt[0] == 'plot_SED':
                try:
                    self.plot_SED = bool(int(txt[txt.index('=') + 1]))
                except:
                    print_error('plot_SED')
                    
            elif txt[0] == 'plot_cutoff_freq':
                try:
                    self.plot_cutoff_freq = int(int(txt[txt.index('=') + 1]))
                except:
                    print_error('plot_cutoff_freq')

            elif txt[0] == 'plot_interval':
                try:
                    self.plot_interval  = int(int(txt[txt.index('=') + 1]))
                except:
                    print_error('plot_interval')

            elif txt[0] == 'plot_slice':
                try:
                    self.plot_slice = bool(int(txt[txt.index('=') + 1]))
                except:
                    print_error('plot_slice')

            elif txt[0] == 'qpoint_slice_index':
                try:
                    self.q_slice_index = int(txt[txt.index('=') + 1])
                except:
                    print_error('qpoint_slice_index')

            ## Lorentz fitting
            elif txt[0] == 'lorentz':
                try:
                    self.lorentz = bool(int(txt[txt.index('=') + 1]))
                except:
                    print_error('lorentz')

            elif txt[0] == 'if_show_figures':
                try:
                    self.if_show_figures = bool(int(txt[txt.index('=') + 1]))
                except:
                    print_error('if_show_figures')

            elif txt[0] == 're_output_total_freq_lifetime':
                try:
                    self.re_output_total_freq_lifetime = bool(int(txt[txt.index('=') + 1]))
                except:
                    print_error('re_output_total_freq_lifetime')

            elif txt[0] == 'peak_height':
                try:
                    self.peak_height = float(txt[txt.index('=') + 1])  # The index of label('=') plus 1
                except:
                    print_error('peak_height')

            elif txt[0] == 'peak_prominence':
                try:
                    self.peak_prominence = float(txt[txt.index('=') + 1])  # The index of label('=') plus 1
                except:
                    print_error('peak_prominence')
                    
            elif txt[0] == 'peak_max_hwhm':
                try:
                    self.peak_max_hwhm = float(txt[txt.index('=') + 1])  # The index of label('=') plus 1
                except:
                    print_error('peak_max_hwhm')

            elif txt[0] == 'lorentz_fit_all_qpoint':
                try:
                    self.lorentz_fit_all_qpoint = bool(int(txt[txt.index('=') + 1]))
                except:
                    print_error('lorentz_fit_all_qpoint')

if __name__ == "__main__":

    my_sample = get_parse_input()

'''
@author:
**************************  LiangTing ***************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/4/26 00:03:21 *********************
'''
import time
import h5py
import numpy as np
import os

def compress(params):

    start = time.perf_counter()
    # check if files exist
    if not os.path.exists(params.vels_file):
        print('\nERROR: file {} not found\n'.format(params.vels_file))
        exit()
    if not os.path.exists(params.pos_file):
        print('\nERROR: file {} not found\n'.format(params.pos_file))
        exit()

    print('\nCompressing velocity and position data into .hdf5 database\n'
          'This may take a while...\n')

    if params.file_format == 'lammps':
        # write to hdf5 database
        num_frame = params.total_num_steps // params.output_data_stride            # number of times actually printed
        # Get the total number of atoms from dump file
        with open(params.vels_file, 'r') as fin:
            for j in range(3):  # Always skip the first three lines
                fin.readline()
            dump_num_atoms = int(fin.readline().strip().split()[0])
        fin.close()

        if dump_num_atoms != params.num_atoms:
            raise ValueError('The number of atoms you input in the control'
                             ' file is not equal to the number of atoms in your dump file !!!')

        with h5py.File(params.output_hdf5, 'w') as fout:
            # For velocity
            with open(params.vels_file, 'r') as fin:
                vels_dataset = fout.create_dataset('velocity', (num_frame, params.num_atoms,
                                                                3))     # Initialize the velocity tensor
                vels = np.zeros((params.num_atoms, 3))
                # look them up in lammps output file
                for i in range(num_frame):
                    for j in range(9):                               # Always skip the first nine lines
                        fin.readline()

                    for k in range(params.num_atoms):
                        vels[k, :] = fin.readline().strip().split()[2:]

                    vels_dataset[i, :, :] = vels
            fin.close()

            # For positions
            with open(params.pos_file, 'r') as fin:
                pos_dataset = fout.create_dataset('position', (num_frame, params.num_atoms,
                                                                3))       # Initialize the velocity tensor
                pos = np.zeros((params.num_atoms, 3))
                # look them up in lammps output file
                for i in range(num_frame):
                    for j in range(9):  # Always skip the first nine lines
                        fin.readline()
                    for k in range(params.num_atoms):
                        pos[k, :] = fin.readline().strip().split()[2:]

                    pos_dataset[i, :, :] = pos
            fin.close()

        fout.close()

    print('Done compressing {} and {} into .hdf5 format.'
         '\nThe compressed file is \'{}\' (DON\'T CHANGE IT!)\n'
         .format(params.vels_file, params.pos_file, params.output_hdf5))

    end = time.perf_counter()

    print ('The time taken to compress the file was ' + str(end-start) + ' s\n')
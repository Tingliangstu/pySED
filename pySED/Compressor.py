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
from sys import exit
import re
import sys

def _read_lammps_block(fin, nblock, num_atoms, data_dtype):
    vels_block = np.zeros((nblock, num_atoms, 3), dtype=data_dtype)
    for bi in range(nblock):
        for _ in range(9):  # Always skip the first nine lines
            fin.readline()
        for k in range(num_atoms):
            vels_block[bi, k, :] = fin.readline().strip().split()[2:]
    return vels_block

def _read_gpumd_block(fin, nblock, num_atoms, data_dtype):
    vels_block = np.zeros((nblock, num_atoms, 3), dtype=data_dtype)
    pos_block = np.zeros((nblock, num_atoms, 3), dtype=data_dtype)
    for bi in range(nblock):
        fin.readline()          # Skip number of atoms line
        fin.readline()          # Skip comment line
        for k in range(num_atoms):
            data = fin.readline().strip().split()
            pos_block[bi, k, :] = data[1:4]
            vels_block[bi, k, :] = data[4:7]
    return pos_block, vels_block

def _print_progress(label, current, total):
    if total <= 0:
        return
    percent = int((current / total) * 100)
    bar_len = 40
    filled = int(bar_len * percent / 100)
    bar = "=" * filled + "-" * (bar_len - filled)
    sys.stdout.write(f"\r{label} [{bar}] {percent}%")
    sys.stdout.flush()
    if current >= total:
        sys.stdout.write("\n")

def compress(params):
    start = time.perf_counter()
    # Now pySED only support lammps or gpumd trajectory, future maybe support AIMD
    num_frame = params.total_num_steps // params.output_data_stride   # number of times actually printed
    use_float32 = bool(getattr(params, "use_float32", 0))             # default is float64
    data_dtype = np.float32 if use_float32 else np.float64
    h5_compression = "lzf"
    block_size = int(getattr(params, "compress_block_size", 0))
    if block_size < 1:
        block_size = max(200, num_frame // 40)
        block_size = min(block_size, 2000)

    box_bounds = []

    if params.file_format == 'lammps':
        # check if files exist
        if not os.path.exists(params.vels_file):
            print('\nERROR: file {} not found\n'.format(params.vels_file))
            exit()
        if not os.path.exists(params.pos_file):
            print('\nERROR: file {} not found\n'.format(params.pos_file))
            exit()

        print('\nCompressing velocity and position data into .hdf5 database\n'
              'This may take a while ...\n')

        # Read box information and number of atoms from the first frame of the position file
        with open(params.pos_file, 'r') as fin:
            # Skip the first three lines (TIMESTEP and headers)
            for _ in range(3):
                fin.readline()

            # The fourth line contains the number of atoms
            dump_num_atoms = int(fin.readline().strip())
            if dump_num_atoms != params.num_atoms:
                raise ValueError('The number of atoms you input in the control file '
                                 'is not equal to the number of atoms in your dump file!')

            # The next line contains BOX BOUNDS header
            box_bounds_header = fin.readline()
            is_tilted = "xy" in box_bounds_header.lower()  # Check if it's a tilted box

            # The next three lines contain box bounds
            for _ in range(3):
                bounds_line = fin.readline().strip().split()
                bounds = [float(x) for x in bounds_line]
                box_bounds.append(bounds)
        fin.close()
        # Rearrange box bounds into matrix format
        if is_tilted:
            # Extract tilted box parameters
            xlo_bound, xhi_bound, xy = box_bounds[0]  # x-dimension boundaries and tilt
            ylo_bound, yhi_bound, xz = box_bounds[1]  # y-dimension boundaries and tilt
            zlo_bound, zhi_bound, yz = box_bounds[2]  # z-dimension boundaries and tilt
            # Adjust xlo and xhi (see https://docs.lammps.org/Howto_triclinic.html)
            xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
            xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
            # Adjust ylo and yhi
            ylo = ylo_bound - min(0.0, yz)
            yhi = yhi_bound - max(0.0, yz)
            # zlo and zhi remain the same
            zlo = zlo_bound
            zhi = zhi_bound
            # Correct calculations for lx, ly, lz
            lx = xhi - xlo  # x-dimension length
            ly = yhi - ylo  # y-dimension length
            lz = zhi - zlo  # z-dimension length
            # Cell matrix for a tilted box
            cell_matrix = [
                [lx, 0.0, 0.0],
                [xy, ly, 0.0],
                [xz, yz, lz]
            ]
        else:
            # Extract orthogonal box parameters
            xlo, xhi = box_bounds[0][:2]  # x-dimension boundaries
            ylo, yhi = box_bounds[1][:2]  # y-dimension boundaries
            zlo, zhi = box_bounds[2][:2]  # z-dimension boundaries

            # Correct calculations for lx, ly, lz
            lx = xhi - xlo  # x-dimension length
            ly = yhi - ylo  # y-dimension length
            lz = zhi - zlo  # z-dimension length

            # Cell matrix for an orthogonal box
            cell_matrix = [
                [lx, 0.0, 0.0],
                [0.0, ly, 0.0],
                [0.0, 0.0, lz]
            ]

        # Convert to numpy array for consistency and further use
        box_bounds = np.array(cell_matrix)

        # get velocity and position information
        with h5py.File(params.output_hdf5, 'w') as fout:
            fout.create_dataset('box', data=box_bounds)  # Store box data
            # For velocity
            with open(params.vels_file, 'r') as fin:
                vels_dataset = fout.create_dataset(
                    'velocity',
                    (num_frame, params.num_atoms, 3),
                    dtype=data_dtype,
                    compression=h5_compression,
                )  # Initialize the velocity tensor
                # look them up in lammps output file
                for i0 in range(0, num_frame, block_size):
                    nblock = min(block_size, num_frame - i0)
                    vels_block = _read_lammps_block(fin, nblock, params.num_atoms, data_dtype)
                    vels_dataset[i0:i0 + nblock, :, :] = vels_block
                    _print_progress("  Compressing lammps velocity", i0 + nblock, num_frame)
            fin.close()

            # For positions
            with open(params.pos_file, 'r') as fin:
                pos_dataset = fout.create_dataset(
                    'position',
                    (num_frame, params.num_atoms, 3),
                    dtype=data_dtype,
                    compression=h5_compression,
                )  # Initialize the velocity tensor
                # look them up in lammps output file
                for i0 in range(0, num_frame, block_size):
                    nblock = min(block_size, num_frame - i0)
                    pos_block = _read_lammps_block(fin, nblock, params.num_atoms, data_dtype)
                    pos_dataset[i0:i0 + nblock, :, :] = pos_block
                    _print_progress("  Compressing lammps position", i0 + nblock, num_frame)
            fin.close()

        fout.close()
        print('\nDone compressing {} and {} into .hdf5 format.'
              '\nThe compressed file is \'{}\' (DON\'T CHANGE IT!)\n'
              .format(params.vels_file, params.pos_file, params.output_hdf5))

    ############################### For xyz ######################################
    if params.file_format == 'gpumd':

        if not os.path.exists(params.dump_xyz_file):
            print('\nERROR: file {} not found\n'.format(params.dump_xyz_file))
            exit()

        print('\nCompressing velocity and position data from GPUMD extxyz file into .hdf5 database\n'
              'This may take a while ...\n')

        with open(params.dump_xyz_file, 'r') as fin:
            first_line = fin.readline()
            num_atoms = int(first_line.strip().split()[0])
            if num_atoms != params.num_atoms:
                raise ValueError(
                    'The number of atoms you input in the control file is'
                    ' not equal to the number of atoms in your extxyz file !!!')
            # Read second line to get box information
            lattice_line = fin.readline().strip()
            if "Lattice=" in lattice_line or "lattice=" in lattice_line:
                # Extract Lattice parameters
                lattice_match = re.search(r'lattice\s*=\s*\"([^\"]+)\"', lattice_line, re.IGNORECASE)
                if not lattice_match:
                    raise ValueError('Lattice parameter not found in the GPUMD file.')
                lattice_data = lattice_match.group(1)
                box_bounds = np.array(lattice_data.split(), dtype=float).reshape(3, 3)
            else:
                raise ValueError('Lattice parameter not found in the GPUMD file.')

        fin.close()

        with h5py.File(params.output_hdf5, 'w') as fout:
            fout.create_dataset('box', data=box_bounds)  # Store box data
            vels_dataset = fout.create_dataset(
                'velocity',
                (num_frame, params.num_atoms, 3),
                dtype=data_dtype,
                compression=h5_compression,
            )
            pos_dataset = fout.create_dataset(
                'position',
                (num_frame, params.num_atoms, 3),
                dtype=data_dtype,
                compression=h5_compression,
            )
            with open(params.dump_xyz_file, 'r') as fin:
                for i0 in range(0, num_frame, block_size):
                    nblock = min(block_size, num_frame - i0)
                    pos_block, vels_block = _read_gpumd_block(fin, nblock, params.num_atoms, data_dtype)
                    pos_dataset[i0:i0 + nblock, :, :] = pos_block
                    vels_dataset[i0:i0 + nblock, :, :] = vels_block
                    _print_progress("  Compressing gpumd trajectory", i0 + nblock, num_frame)

            fin.close()
        fout.close()
        
        print('\nDone compressing {} into .hdf5 format.'
              '\nThe compressed file is \'{}\' (DON\'T CHANGE IT!)\n'
              .format(params.dump_xyz_file, params.output_hdf5))

    end = time.perf_counter()
    print('The time taken to compress the file was ' + str(end - start) + ' s\n')

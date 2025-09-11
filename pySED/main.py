#!/usr/bin/env python

'''
@author:
**************************  LiangTing ***************************
                    liangting.zj@gmail.com
************************ 2024/12/07 23:03:21 ********************
'''

# Python modules
import sys
import os
import h5py

# My modules
from pySED import (My_Parsers,
                   Compressor,
                   construct_BZ,
                   Phonon,
                   FileIO,
                   Plot_SED,
                   Lorentz)
                   
def main():
    # ******************************** Show pySED imformation ************************
    logo_lines = [
     "                    _ __   _   _    _____  ______  _____",
     "                   | '_ \\ | | | |  / ____||  ____||  __ \\",
     "                   | |_) || |_| | | (___  | |__   | |  | |",
     "                   | .__/  \__, |  \___ \ |  __|  | |  | |",
     "                   | |      __/ | ____) | | |____ | |__| |",
     "                   |_|     |___/ |_____/  |______||_____/",
     "                                              v-2.2.0"
     ]
    for line in logo_lines:
        print(line)
    print('***************** Now, pySED can interface with GPUMD/LAMMPS **************************')
    print('******* To implement the SED method in 2010, phonon lifetime can be calculated ********')
    print('************************ Author: liangting.zj@gmail.com *******************************')
    print('***************************************************************************************')

    # ********************************  Get the input file name ********************************
    if len(sys.argv) == 1:
        input_file = 'input_SED.in'
    elif len(sys.argv) > 2:
        print('\nERROR: pySED takes 1 or 0 arguments!'
                '\nTry \'pySED --help or pySED -h\' for more information.\n')
        exit()

    # ******************************** For usage help ********************************
    elif len(sys.argv) == 2 and sys.argv[1].lower() in ['h', 'help', '-h', '--help', '-H', '-Help']:
        help_text = r"""
    USAGE:
        pySED [input_file] or pysed or pySED

    DESCRIPTION:
        pySED is a tool for calculating Phonon Spectral Energy Density (SED) 
        and phonon lifetime from Molecular Dynamics (MD) simulation trajectories.
        It supports GPUMD/LAMMPS trajectories output formats.

        If no input_file is given, defaults to `input_SED.in` in current directory.
        
    REFERENCE:
    [1] Journal of Applied Physics 138, 075101 (2025).
        https://doi.org/10.1063/5.0215411
        
    [2] Phys. Rev. B 81, 081411 (2010).
        https://doi.org/10.1103/PhysRevB.81.081411
        
    Please cite above papers if you use pySED package in your research.

    INPUT PARAMETERS (input_SED.in):

    ────────────────────────────────────────────────────────────
    [ MD Simulation Parameters ]
        num_atoms               : Number of atoms in simulation structure. 
        total_num_steps         : MD simulation total steps.
        time_step               : Time step in femtoseconds (fs). 
        output_data_stride      : Data output stride in steps.

    ────────────────────────────────────────────────────────────
    [ Data Splitting & Compression ]
        num_splits              : Number of blocks to average. Default = 1. Tips: set num_splits = 5~10 is better.
        compress                : Compress velocity/position data to HDF5. Default = 1. Tips: do not change it.

    ────────────────────────────────────────────────────────────
    [ Input / Output Files ]
        out_files_name          : Prefix for output files. Customize your output file name.
        basis_lattice_file      : Basis lattice file path. Default = 'basis.in'.
        file_format             : 'gpumd' or 'lammps'. Default = 'gpumd'. Tips: set it to 'lammps' if you use LAMMPS.
        
        (For GPUMD user) 
        dump_xyz_file           : Dump file GPUMD's trajectories format. Default = 'dump.xyz'.
        
        (For LAMMPS user) 
        pos_file                : Atomic positions file for LAMMPS. Default = 'pos.dat'.
        vels_file               : Atomic velocities file for LAMMPS. Default = 'vels.dat'.
        lammps_unit             : LAMMPS unit style ('metal' or 'real'). Default = 'metal'.
        
        output_hdf5             : Compressed HDF5 file name. Default = 'vel_pos_compress.hdf5'.

    ────────────────────────────────────────────────────────────
    [ Parallelization ]
        use_parallel            : Enable multithreading (1=yes, 0=no). Default = 1.
        max_cores               : Max CPU cores for parallel mode. Default = 4. Tips: set max_cores = 4, 6, 8 is better.

    ────────────────────────────────────────────────────────────
    [ Crystal Structure ]
        prim_unitcell           : Primitive unit cell parameters, one can check it use ovito (POSCAR).
        supercell_dim           : Supercell dimensions. Default = [1,1,1]. Tips: the larger the supercell, the denser the commensurate q-point sampling.
        rescale_prim            : Rescale primitive cell (1=yes, 0=no). Default = 1. Tips: for structural optimization using the NPT ensemble.

    ────────────────────────────────────────────────────────────
    [ Q-points & BZ ]
        num_qpaths              : Number of dispersion paths. Default = 1
        q_path_name             : Name of q-path. Default='GA'. Tips: for gaphene, one can set q_path_name = 'GMKG'.
        q_path                  : Q-points in fractional coordinates. Tips: for q_path = 'GA', q_path = 0.0 0.0 0.0  0.5 0.0 0.0.

    ────────────────────────────────────────────────────────────
    [ Plotting Options ]
        plot_SED                : Plot/fitting mode (1) or compute mode (0). Default: plot_SED = 0.
        plot_cutoff_freq        : Max frequency for plot (THz). Default = None. Tips: if None, the max frequency is the max frequency in the SED data.
        plot_interval           : y-axis tick interval (THz). Default = 5.
        qpoint_slice_index      : Index of q-point for single plot. Tips: start from 0.
        plot_slice              : Plot q-slice (1=yes, 0=no). Default = 0.
        if_show_figures         : Show figures on screen. Default = 0.

    ────────────────────────────────────────────────────────────
    [ Lorentzian Fitting ] (We highly encourage users to read through the tips below carefully)
    
        lorentz                 : Enable Lorentzian fitting (1=yes, 0=no). Default = 0.
        peak_height             : Minimum peak height. Default = None.
        peak_prominence         : Minimum peak prominence. Default = None.
        initial_guess_hwhm      : Initial guess HWHM. Default = 0.001.
        lorentz_fit_cutoff      : Frequency cutoff for fitting (THz). Default = None.
        
        lorentz_fit_all_qpoint  : Fit all q-points. Default = 0. Tips: fit all q-points peak and output the phonon lifetime (ps).
        re_output_total_freq_lifetime : Recalculate lifetimes for typical q-point by set `qpoint_slice_index`. Default = 0.
        
        
        Tips for lorentzian fitting：
        1. `peak_height`: the minimum height of a peak, which is the minimum height threshold of an SED peak to be taken into account.
           If `peak_height` is too large, SED peaks may be missed; if too small, noise might be counted as peaks.
           You can set `qpoint_slice_index = xx` (set q-point index xx) and `plot_slice = 1`
           to determine the SED peak height and prominence.
           
        2. The prominence of a peak measures how much a peak stands out from the surrounding baseline of the signal,
           and is defined as the vertical distance between the peak and its lowest contour line.
           
        3. If you find that the peak fitting for certain q-points is unsatisfactory, 
           you can set `qpoint_slice_index = xx` and `re_output_total_freq_lifetime = 1` to readjust them.
           In this case, be sure to set `lorentz_fit_all_qpoint = 0`.
           
           You can adjust `peak_height` and `peak_prominence`, 
           and pySED will reoutput to the `TOTAL-LORENTZ-Qpoints.Fre_lifetime` file.

        
    ────────────────────────────────────────────────────────────
    EXAMPLES:
        Run SED calculation:
            pySED input_SED.in

        Plot/Fitting existing SED:
            pySED input_SED.in     # with plot_SED = 1 (plot) / plot_SED = 1 and lorentz = 1 (Fitting)

        Show help:
            pySED -h or pysed -h
    """
        print(help_text)
        exit()
    else:
        input_file = str(sys.argv[1])

    # ******************************** Check if file exist ********************************
    if not os.path.exists(input_file):
        print('\n**************** ERROR: File \'{}\' not found in current path! ****************'.format(input_file))
        print('\n*********************** Please check the file name and path. **************************')
        print('\n******* You can run "pysed -h" to see usage and example `input_SED.in` format. ********\n')
        exit()

    ### **************** Read the input parameter file and get the control parameters ******************
    params = My_Parsers.get_parse_input(input_file)

    ## **************** Select the SED-modes (compute or plot-including-fitting mode) ******************
    if params.plot_SED:
        pass

    else:
        ### *********************** Compress the velocity and position data to hdf5  ***********************
        print('\n******************* You are in the SED-Computing mode ********************')
        if params.compress:
            # print(params.output_hdf5 + ' file already exists, and there is no need to compression !!')
            if os.path.exists(params.output_hdf5):
                print('\n*************** ' + params.output_hdf5 + ' file will be used ! ****************')
            else:
                Compressor.compress(params)
                print('******** Compress the velocity and position data to hdf5 is DONE! ********')
        else:
            raise EOFError('**** Please set "compress = 1" in input file to compress the velocity and position data ****')

        # ***************************** For implement SED Method *******************************
        ## Construct BZ
        BZ_lattice_info = construct_BZ.BZ_methods(params)

        ## compute SED
        phonons = Phonon.spectral_energy_density(params)
        phonons.compute_sed(params, BZ_lattice_info)

        ### Save the data to output files
        FileIO.write_output(phonons, params, BZ_lattice_info)

        print('\n******* SED-Computing mode is running done for all Q-points, now you can run plotting mode *******\n')
        print('******************** Please set "plot_SED = 1" in input file to plot the SED. ********************\n')
        exit()

    # ***************************** Finished SED compute, For plotting *******************************
    # ********** If found params.lorentz_fit_all_qpoint, just fitting all Q points *******************

    if params.lorentz_fit_all_qpoint and params.lorentz:
        data = FileIO.load_data(params)
        print('\n********** You are in the ALL plotting mode (total {} Qpoint to fitting) !!! **********\n'
              .format(len(data.q_distances)))
        for j in range(len(data.q_distances)):
            # Control parameters for all Qpoints
            params.q_slice_index = j
            params.if_show_figures = False
            params.plot_lorentz = False

            print('\n***************************************************************************************')
            print('***************** Now processing {}-th Qpoint SED Lorentz fitting **********************'.format(j))
            data = FileIO.load_data(params)
            #Plot_SED.plot_slice(data, params)              # no need to output here
            Lorentz.lorentz(data, params)
            print('***************************************************************************************')

        FileIO.deal_total_fre_lifetime(params, len(data.q_distances))
        print('\n****************************** pySED PLOTTING ALL DONE !!!!! ****************************')
        print('*****************************************************************************************')
        print(' If some Q-points don\'t fit well enough, one can rerun the single fitting mode by setting:')
        print(' \"qpoint_slice_index = XX, re_output_total_freq_lifetime = 1, lorentz_fit_all_qpoint = 0\"')
        print('*****************************************************************************************')
        exit()

    # ************************* Doing lorentz fitting for single Q-point *****************************
    if params.plot_SED:
        if params.plot_SED:
            print('\n***** You are in the single SED-qpoint plotting mode, plot SED for single Q-point *****\n')
            data = FileIO.load_data(params)
            Plot_SED.plot_bands(data, params)

            if params.plot_slice:
                Plot_SED.plot_slice(data, params)

         ### Fit to a Lorentz function
        if params.lorentz:
           print('\n******************* You are doing lorentz fitting for single Q-point ******************\n')
           data = FileIO.load_data(params)
           Lorentz.lorentz(data, params)
           print('\n************************** Single Lorentz fitting ALL DONE! ***************************')

        # This parameter allows you to re-adjust the single peak of a single Q-point and re-output the total phonon lifetime
        if params.re_output_total_freq_lifetime:
            data = FileIO.load_data(params)
            FileIO.deal_total_fre_lifetime(params, len(data.q_distances))

        print('\n************************** pySED SINGLE-PLOTTING ALL DONE !!!!! ***********************\n')
        print('***************************************************************************************')
        print('  If single Q-point fitting works well, one can fit all Q-point directly by setting:')
        print('                              \"lorentz_fit_all_qpoint = 1\"')
        print('***************************************************************************************\n')

if __name__ == '__main__':
    main()

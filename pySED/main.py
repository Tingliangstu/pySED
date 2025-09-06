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
     "                                              v-2.1.0"
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
    elif len(sys.argv) == 2 and (sys.argv[1] == 'h' or sys.argv[1] == 'help' or
            sys.argv[1] == 'HELP' or sys.argv[1] == '-h' or sys.argv[1] == '--help' or
            sys.argv[1] == '--HELP'):

        print('\nUSAGE: pySED [input_file]\n\n\'input_file\' should be the name of the file '
            'containing the parameters \nto calculate the Phonon Spectral Energy Density (SED).'
            '\n\nIf no input_file name is given, the default is \'input_SED.in\'.\n')
        exit()
    else:
        input_file = str(sys.argv[1])

    # ******************************** Check if file exist ********************************
    if not os.path.exists(input_file):
        print('\n********************* ERROR: File \'{}\' not found in current path! *********************\n'.format(input_file))
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

        print('\n********** SED-Computing mode is running done for all Q-points, Now you can run plotting mode **********\n')
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

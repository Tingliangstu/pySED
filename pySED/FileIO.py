'''
@author:
**************************  LiangTing ****************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/5/15 23:03:21 **********************
'''

import os
import numpy as np

def write_output(phonons, params, BZ_lattice_info):
    output_partial = getattr(params, 'output_partial', 0)
    if output_partial:
        # Overall SED = sum over atom type and direction
        sed_total = phonons.sed_avg.sum(axis=(-2, -1))
        np.savetxt(params.out_files_name + '.SED', sed_total)

        partial_dir = params.out_files_name + '_partial_SED'
        os.makedirs(partial_dir, exist_ok=True)
        dir_labels = ['x', 'y', 'z']

        for t_idx in range(phonons.sed_avg.shape[2]):
            for d_idx, d_lab in enumerate(dir_labels):
                out_path = os.path.join(partial_dir, f"{params.out_files_name}.SED_type{t_idx+1}_{d_lab}")
                np.savetxt(out_path, phonons.sed_avg[:, :, t_idx, d_idx])
    else:
        np.savetxt(params.out_files_name + '.SED', phonons.sed_avg)

    np.savetxt(params.out_files_name + '.Qpts', BZ_lattice_info.reduced_qpoints, fmt='%.8f')
    np.savetxt(params.out_files_name + '.THz', phonons.freq_fft, fmt='%.8f')

    with open(params.out_files_name + '.Q_distances_and_labels', 'w') as f:
        # q_distances
        f.write("Global distances along the paths:\n")
        f.write(" ".join(f"{d:.10f}" for d in BZ_lattice_info.q_distances) + "\n\n")

        # q_labels
        f.write("High-symmetry points and their distances:\n")
        for distance, label in BZ_lattice_info.q_labels:
            f.write(f"{float(distance):.10f}   {label}\n")

class load_data(object):

    def __init__(self, params):

        if getattr(params, 'plot_partial_SED', 0):

            partial_dir = params.out_files_name + '_partial_SED'
            type_idx = params.plot_partial_type

            if type_idx is None or type_idx < 0:
                raise ValueError('\n*************** plot_partial_SED type index is invalid ***************')

            if params.plot_partial_dir is None:
                # Sum x/y/z for the given type
                file_x = os.path.join(partial_dir, f"{params.out_files_name}.SED_type{type_idx+1}_x")
                file_y = os.path.join(partial_dir, f"{params.out_files_name}.SED_type{type_idx+1}_y")
                file_z = os.path.join(partial_dir, f"{params.out_files_name}.SED_type{type_idx+1}_z")
                if not (os.path.exists(file_x) and os.path.exists(file_y) and os.path.exists(file_z)):
                    raise FileNotFoundError(
                        '\n*************** plot_partial_SED type index is out of range or files missing ***************')
                sed_x = np.loadtxt(file_x)
                sed_y = np.loadtxt(file_y)
                sed_z = np.loadtxt(file_z)
                self.sed_avg = sed_x + sed_y + sed_z
            else:
                d = params.plot_partial_dir
                file_d = os.path.join(partial_dir, f"{params.out_files_name}.SED_type{type_idx+1}_{d}")
                if not os.path.exists(file_d):
                    raise FileNotFoundError(
                        '\n*************** plot_partial_SED type index is out of range or file missing ***************')
                self.sed_avg = np.loadtxt(file_d)

        else:
            self.sed_avg = np.loadtxt(params.out_files_name + '.SED')

        self.qpoints = np.loadtxt(params.out_files_name + '.Qpts')
        self.freq_fft = np.loadtxt(params.out_files_name + '.THz')

        self.q_distances = []
        self.q_labels = {}
        with open(params.out_files_name + '.Q_distances_and_labels', 'r') as f:
            lines = f.readlines()
            # q_distances
            if lines[0].strip() == "Global distances along the paths:":
                self.q_distances = [float(x) for x in lines[1].split()]
            # q_labels
            if lines[3].strip() == "High-symmetry points and their distances:":
                for line in lines[4:]:
                    distance, label = line.split(maxsplit=1)
                    self.q_labels[float(distance)] = label.strip()

def write_lorentz(lorentz, params):
    np.savetxt(params.out_files_name + '_LORENTZ-{}.params'.format(params.q_slice_index), lorentz.popt)
    np.savetxt(params.out_files_name + '_LORENTZ-{}.error'.format(params.q_slice_index), lorentz.pcov)
    print('**************** The specific parameters of the fit are successfully written ***************')

def write_phonon_lifetime(lorentz, params):

    out_lifetime_file = 'Generate by pySED codes, Email: liangting.zj@gmail.com\n'
    out_lifetime_file += "First_line: Frequency (THz)  Second_line: Phonon Lifetime (ps)\n"

    for i in range(len(lorentz.popt)):

        if lorentz.popt[i][2] == 0:  # don't output fitting fail frequency
            continue

        out_lifetime_file += '{0:.6f} {1:.8f} \n'.format(lorentz.popt[i][0], 1/(2 * np.pi * lorentz.popt[i][2]))
        # write the file

    f = open('LORENTZ-{}-th-Qpoints.Fre_lifetime'.format(params.q_slice_index), 'w')
    f.write(out_lifetime_file)
    f.close()

    print('************** LORENTZ-{}-th Qpoints.Fre_lifetime is written successfully **************'.format(params.q_slice_index))

def deal_total_fre_lifetime(params, total_qpoints):

    out_lifetime_file = 'Generate by pySED codes, Email: liangting.zj@gmail.com\n'
    out_lifetime_file += "First_line: Frequency (THz)  Second_line: Phonon Lifetime (ps)\n"

    total_num_Fre_lifetime = 0

    for i in range(total_qpoints):
        load_file_name = 'LORENTZ-{}-th-Qpoints.Fre_lifetime'.format(i)
        try:
            with open(load_file_name, 'r') as f:
                lines = [line.strip() for line in f.readlines()]
                # Check if the file only contains the header lines
                if len(lines) <= 2 or all(line == '' for line in lines[2:]):
                    print(f'\nWarning: {load_file_name} does not contain fitted phonon lifetimes, skipping.')
                    continue

            # Load the numerical data
            Freq, lifetime = np.loadtxt(load_file_name, skiprows=2, unpack=True)

            if isinstance(Freq, np.float64):
                Freq = np.array([Freq])
                lifetime = np.array([lifetime])

            for j in range(len(Freq)):
                out_lifetime_file += '{0:.6f} {1:.8f}\n'.format(Freq[j], lifetime[j])
                total_num_Fre_lifetime += 1

        except:
            raise FileNotFoundError(
                '\n*************** File LORENTZ-{}-th-Qpoints.Fre_lifetime reading ERROR ***************'.format(i))

    f = open('TOTAL-LORENTZ-Qpoints.Fre_lifetime', 'w')
    f.write(out_lifetime_file)
    f.close()

    if not params.re_output_total_freq_lifetime:
        print('\n**** TOTAL-LORENTZ-Qpoints.Fre_lifetime is written successfully (Total {} points) *****'
                                                                                    .format(total_num_Fre_lifetime))
    else:
        print('\n*** TOTAL-LORENTZ-Qpoints.Fre_lifetime is Re-written successfully (Total {} points) ***'
              .format(total_num_Fre_lifetime))

'''
@author:
**************************  LiangTing ****************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/5/15 23:03:21 **********************
'''

import numpy as np

def write_output(phonons, params, lattice):
    np.savetxt(params.out_files_name + '.SED', phonons.sed_avg)
    np.savetxt(params.out_files_name + '.Qpts', lattice.reduced_qpoints, fmt='%.8f')
    np.savetxt(params.out_files_name + '.THz', phonons.freq_fft, fmt='%.8f')

    with open(params.out_files_name + '.Q_distances_and_labels', 'w') as f:
        # q_distances
        f.write("Global distances along the paths:\n")
        f.write(" ".join(f"{d:.8f}" for d in lattice.q_distances) + "\n\n")

        # q_labels
        f.write("High-symmetry points and their distances:\n")
        for distance, label in lattice.q_labels.items():
            f.write(f"{float(distance):.8f} {label}\n")

class load_data(object):

    def __init__(self, params):

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

        out_lifetime_file += '{0:.6f} {1:.8f} \n'.format(lorentz.popt[i][0], 1/(2 * lorentz.popt[i][2]))
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
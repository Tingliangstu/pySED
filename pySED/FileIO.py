'''
@author:
**************************  LiangTing ****************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/5/15 23:03:21 **********************
'''

import numpy as np

def write_output(phonons, params, lattice):
    np.savetxt(params.out_files_name + '.SED', phonons.sed_avg)
    np.savetxt(params.out_files_name + '.Qpts', lattice.reduced_qpoints, fmt='%.6f')
    np.savetxt(params.out_files_name + '.THz', phonons.freq_fft, fmt='%.6f')

class load_data(object):

    def __init__(self, params):

        self.sed_avg = np.loadtxt(params.out_files_name + '.SED')
        self.qpoints = np.loadtxt(params.out_files_name + '.Qpts')
        self.freq_fft = np.loadtxt(params.out_files_name + '.THz')

def write_lorentz(lorentz, params):
    np.savetxt(params.out_files_name + '_LORENTZ-{}.params'.format(params.q_slice_index), lorentz.popt)
    np.savetxt(params.out_files_name + '_LORENTZ-{}.error'.format(params.q_slice_index), lorentz.pcov)
    print('**************** The specific parameters of the fit are successfully written ***************')

def write_phonon_lifetime(lorentz, params):

    out_lifetime_file = 'Generate by LT\'s codes, Email: liangting.zj@gmail.com\n'
    out_lifetime_file += "First_line: Frequency (THz)  Second_line: Phonon Lifetime (ps)\n"

    for i in range(len(lorentz.popt)):

        if lorentz.popt[i][2] == 0:  # don't output fitting fail frequency
            continue

        out_lifetime_file += '{0:.4f} {1:.8f} \n'.format(lorentz.popt[i][0], 1/(2 * lorentz.popt[i][2]))
        # write the file

    f = open('LORENTZ-{}-th-Qpoints.Fre_lifetime'.format(params.q_slice_index), 'w')
    f.write(out_lifetime_file)
    f.close()

    print('************** LORENTZ-{}-th Qpoints.Fre_lifetime is written successfully **************'.format(params.q_slice_index))

def deal_total_fre_lifetime(params):

    out_lifetime_file = 'Generate by LT\'s codes, Email: liangting.zj@gmail.com\n'
    out_lifetime_file += "First_line: Frequency (THz)  Second_line: Phonon Lifetime (ps)\n"

    total_num_Fre_lifetime = 0
    for i in range(sum(params.num_qpoints)):

        try:
            load_file_name = 'LORENTZ-{}-th-Qpoints.Fre_lifetime'.format(i)
            Freq, lifetime = np.loadtxt(load_file_name, skiprows=2, unpack=True)

            for j in range(len(Freq)):
                out_lifetime_file += '{0:.4f} {1:.8f}\n'.format(Freq[j], lifetime[j])
                total_num_Fre_lifetime += 1

        except:
            raise FileNotFoundError('*************** File LORENTZ-{}-th-Qpoints.Fre_lifetime reading ERROR ***************'
                                    .format(i))

    f = open('TOTAL-LORENTZ-Qpoints.Fre_lifetime', 'w')
    f.write(out_lifetime_file)
    f.close()

    if not params.re_output_total_freq_lifetime:
        print('**** TOTAL-LORENTZ-Qpoints.Fre_lifetime is written successfully (Total {} points) *****'
                                                                                    .format(total_num_Fre_lifetime))
    else:
        print('\n*** TOTAL-LORENTZ-Qpoints.Fre_lifetime is Re-written successfully (Total {} points) ***'
              .format(total_num_Fre_lifetime))
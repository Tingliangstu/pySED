'''
@author:
**************************  LiangTing ****************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/5/15 23:03:21 **********************
'''

import numpy as np
from pySED import FileIO, Plot_SED
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

class lorentz:
    def __init__(self, data, params):

        def lorentzian(xarr, center, amplitude, hwhm):
            return  amplitude / (1 + ((xarr - center) / hwhm) ** 2)

        self.q_index = params.q_slice_index
        self.lorentz_fit_cutoff = params.lorentz_fit_cutoff
        self.sed = data.sed_avg[:, self.q_index]
        self.thz = data.freq_fft

        if self.lorentz_fit_cutoff:
            index = 0
            for i, x in enumerate(self.thz):
                if x > self.lorentz_fit_cutoff:
                    index = i
                    break
            self.sed = self.sed[0:index]
            self.thz = self.thz[0:index]
            print("****************** WARNING: Now the Frequency will cutoff to {0:5.3f} ********************\n".
                                                                                         format(max(self.thz)))

        ## ********* Using the scipy.signal.find_peaks find the peak in SED curve *********
        peaks, amp_max = find_peaks(self.sed, height=params.peak_height, prominence=params.peak_prominence)

        print('*** Found {} peaks in the SED-{}-th qpoint curve, Please compare with the actual peak ***'.format(len(peaks),
                                                                                                        params.q_slice_index))

        print('**** Peaks is as follows, one can tune fitting paras (peak_height) according them ****:\n\n {} \n'.format(self.sed[peaks]))

        # some bounds on the fitting. Might need to tweak these
        dx = 1                # The size of the peak left and right offset during fitting
        maxfev = 1e15         # Should be the maximum number of fits

        self.xarr = np.arange(len(self.sed))
        self.popt = np.zeros((len(peaks), 3))
        self.pcov = np.zeros((len(peaks), 3))

        #params.bounds = np.zeros((len(peaks), 2))

        for i in range(len(peaks)):
            # Peak started index  # Record boundary
            adjust_number = params.modulate_factor   # for better fitting (next version)
            start = amp_max['left_bases'][i] + adjust_number
            end = amp_max['right_bases'][i] - adjust_number

            #start = peaks[i]-5
            #end = peaks[i]+5            # For dubug
            # Boundary for three parmas
            lb = [self.thz[peaks[i] - dx], amp_max['peak_heights'][i], 1e-12]
            ub = [self.thz[peaks[i] + dx], amp_max['peak_heights'][i]*2, params.peak_max_hwhm]  # np.inf

            ''' p0 = Initial guess for the parameters (length N)
                popt = Optimal values for the parameters
                pcov = Standard deviation errors on the parameters
            '''
            try:
                self.popt[i, :], pcov = curve_fit(lorentzian,
                                                  self.thz[start:end],
                                                  self.sed[start:end],
                                                  p0 = [self.thz[peaks[i]], self.sed[peaks[i]], 1],
                                                  bounds=(lb, ub), maxfev=maxfev)

                self.pcov[i, :] = np.sqrt(np.diag(pcov))    # Get the standard deviation errors

            except:
                print('\nWARNING: Lorentz fit for Frequency-{} failed, maybe due to your unreasonable fitting para setting\n'
                      .format(self.thz[peaks[i]]))
                continue
        # Write to the files
        FileIO.write_phonon_lifetime(self, params)

        # Control params
        params.popt = self.popt
        params.plot_lorentz = True
        Plot_SED.plot_slice(data, params)



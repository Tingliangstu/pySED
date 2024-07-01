'''
@author:
**************************  LiangTing ****************************
        liangting.zj@gmail.com --- Refer from Ty Sterling's script
************************ 2021/5/15 23:03:21 **********************
'''

import numpy as np
from pySED import FileIO, Plot_SED
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, windows
import matplotlib.pyplot as plt

class lorentz:
    def __init__(self, data, params):

        """
        Here for test params:
        params.apply_smoothing = False
        params.modulate_factor and test
        """

        def lorentzian(xarr, center, amplitude, hwhm):
            return amplitude / (1 + ((xarr - center) / hwhm) ** 2)

        def smooth_data(y, window_len=3):  # (hanning windows)
            window = windows.hann(window_len)
            y_padded = np.pad(y, (window_len // 2, window_len // 2), mode='reflect')
            y_smooth = np.convolve(y_padded, window / window.sum(), mode='valid')
            return y_smooth

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

        # smoothing data (hanning windows)
        params.apply_smoothing = False
        if params.apply_smoothing:
            self.sed = smooth_data(self.sed)
            print("****************** WARNING: using smoothing SED data for fitting. ********************\n")

        ## ********* Using the scipy.signal.find_peaks find the peak in SED curve *********
        peaks, amp_max = find_peaks(self.sed, height=params.peak_height, prominence=params.peak_prominence)

        print('*** Found {} peaks in the SED-{}-th qpoint curve, Please compare with the actual peak ***'.format(
            len(peaks),
            params.q_slice_index))

        print(
            '**** Peaks is as follows, one can tune fitting paras according them ****:\nFrequency: {0} THz\nPeak_height: {1} J*s\n'.format(
                self.thz[peaks], self.sed[peaks]))

        # some bounds on the fitting. Might need to tweak these
        dx = 1  # The size of the peak left and right offset during fitting
        maxfev = 1e15  # Should be the maximum number of fits

        self.xarr = np.arange(len(self.sed))
        self.popt = np.zeros((len(peaks), 3))
        self.pcov = np.zeros((len(peaks), 3))

        # params.bounds = np.zeros((len(peaks), 2))

        for i in range(len(peaks)):
            # Peak started index  # Record boundary
            adjust_number = params.modulate_factor  # for better fitting (next version)
            test = True   # for test
            if params.modulate_factor > 0 and test:
                print("Shrinking the fitting range of the X-axis toward the middle for {0}-th peak,"
                      " from [{1:.2f} {2:.2f}] to [{3:.2f} {4:.2f}] THz".format(i,
                                                                         self.thz[amp_max['left_bases'][i]],
                                                                         self.thz[amp_max['right_bases'][i]],
                                                                         self.thz[amp_max['left_bases'][i]+adjust_number],
                                                                         self.thz[amp_max['right_bases'][i]+adjust_number]))

            start = amp_max['left_bases'][i] + adjust_number
            end = amp_max['right_bases'][i] - adjust_number

            # start = peaks[i]-5
            # end = peaks[i]+5              # For dubug
            # Boundary for three parmas
            lb = [self.thz[peaks[i] - dx], amp_max['peak_heights'][i], 1e-14]
            ub = [self.thz[peaks[i] + dx], amp_max['peak_heights'][i] * 2, params.peak_max_hwhm]   # np.inf

            # Initial guesses for fitting parameters
            p0 = [self.thz[peaks[i]], self.sed[peaks[i]], params.initial_guess_hwhm]  # Adjusted initial guess for HWHM

            ''' p0 = Initial guess for the parameters (length N)
                popt = Optimal values for the parameters
                pcov = Standard deviation errors on the parameters
            '''
            try:
                self.popt[i, :], pcov = curve_fit(lorentzian,
                                                  self.thz[start:end],
                                                  self.sed[start:end],
                                                  p0=p0,
                                                  bounds=(lb, ub), maxfev=maxfev)

                self.pcov[i, :] = np.sqrt(np.diag(pcov))  # Get the standard deviation errors

            except:
                print(
                    '\nWARNING: Lorentz fit for Frequency: {} failed, maybe due to your unreasonable fitting para setting\n'
                    .format(self.thz[peaks[i]]))
                continue

        # Write to the files
        FileIO.write_phonon_lifetime(self, params)

        # Control params
        params.popt = self.popt

        print("\n**** The fitting hwhm is {}, Please check it. ****".format(self.popt[:, 2]))
        
        params.plot_lorentz = True
        Plot_SED.plot_slice(data, params)

        ########################## For test smooth (fro test)
        # Plot the data and the fits
        #plt.plot(self.thz, aa, label='Original Data')
        #if params.apply_smoothing:
        #    plt.plot(self.thz, self.sed, label='Smoothed Data')

        #for popt in self.popt:
        #    plt.plot(self.thz, lorentzian(self.thz, *popt), label=f'Fit: center={popt[0]::.3f}')

        #plt.xlabel('Frequency (THz)')
        #plt.ylabel('log(Φ(ω))')
        #plt.legend()
        #plt.show()

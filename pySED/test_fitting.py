import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, windows
from scipy.special import wofz
import matplotlib.pyplot as plt

# Voigt function definition
def voigt(x, center, amplitude, sigma, gamma):
    z = ((x - center) + 1j*gamma) / (sigma * np.sqrt(2))
    return amplitude * np.real(wofz(z)) / (sigma * np.sqrt(2*np.pi))

# Hanning window smoothing function
def smooth_data(y, window_len=5):
    window = windows.hann(window_len)
    y_padded = np.pad(y, (window_len // 2, window_len // 2), mode='reflect')
    y_smooth = np.convolve(y_padded, window / window.sum(), mode='valid')
    return y_smooth

# Generate synthetic data
x_data = np.linspace(0, 2, 400)
y_data = (
    voigt(x_data, 0.75, 1e-7, 0.05, 0.05) +
    voigt(x_data, 1.25, 1e-7, 0.05, 0.05) +
    1e-11 * np.random.normal(size=x_data.size)
)

# Option to smooth data
apply_smoothing = False
if apply_smoothing:
    y_data = smooth_data(y_data)

# Find peaks in the data
peaks, _ = find_peaks(y_data, height=1e-10, prominence=1e-10)

# Prepare for fitting
popt_list = []
pcov_list = []
fit_range = 20  # Adjust this value based on your data

for peak in peaks[:2]:  # Limit to the first two peaks for demonstration
    # Define fitting range around the peak
    start = max(0, peak - fit_range)
    end = min(len(x_data), peak + fit_range)

    # Initial guesses for fitting parameters
    p0 = [x_data[peak], y_data[peak], 0.05, 0.05]  # Initial guess for sigma and gamma

    # Bounds for fitting parameters
    bounds = ([x_data[start], 0, 0, 0], [x_data[end], np.inf, np.inf, np.inf])

    # Perform curve fitting
    try:
        popt, pcov = curve_fit(voigt, x_data[start:end], y_data[start:end], p0=p0, bounds=bounds)
        popt_list.append(popt)
        pcov_list.append(pcov)
    except RuntimeError:
        print(f"Fit for peak at {x_data[peak]} failed.")

# Plot the data and the fits
plt.plot(x_data, y_data, label='Smoothed Data', color='gray')

for popt in popt_list:
    plt.plot(x_data, voigt(x_data, *popt), label=f'Fit: center={popt[0]:.3f}', linewidth=2)

plt.xlabel('Frequency (THz)')
plt.ylabel('log(a(w))')
plt.legend()
plt.title('Voigt Fit')
plt.show()

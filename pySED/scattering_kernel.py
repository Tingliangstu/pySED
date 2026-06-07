"""Shared low-level kernels for pySED scattering calculations.

The functions here are deliberately small CPU/NumPy kernels.  They centralize
the operations that will later be natural targets for optional CuPy or
Numba-CUDA implementations: phase construction, density amplitudes, time FFT
power spectra, and block statistics.
"""

from dataclasses import dataclass

import numpy as np
from scipy.fft import fft, fftfreq


@dataclass
class KernelBackend:
    name: str
    xp: object
    fft: object
    fftfreq: object
    asnumpy: object


@dataclass
class FrequencyGrid:
    frequencies_thz: np.ndarray
    freq_slice: slice


def resolve_backend(backend="cpu"):
    """Resolve a scattering kernel backend.

    ``"cpu"`` and ``"numpy"`` use NumPy/SciPy and are always available.
    ``"cupy"`` and ``"gpu"`` request CuPy/CuPyX.  CuPy is optional and must be
    installed by the user with the wheel that matches their CUDA runtime.
    """

    if isinstance(backend, KernelBackend):
        return backend
    if backend is None or backend in ("cpu", "numpy"):
        return KernelBackend(
            name="cpu",
            xp=np,
            fft=fft,
            fftfreq=fftfreq,
            asnumpy=lambda array: np.asarray(array),
        )
    if backend in ("cupy", "gpu"):
        try:
            import cupy as cp
            from cupyx.scipy.fft import fft as cupy_fft
            from cupyx.scipy.fft import fftfreq as cupy_fftfreq
        except ImportError as exc:
            raise ImportError(
                "CuPy backend requested but CuPy is not installed. Install the "
                "CuPy package matching your CUDA runtime, or use backend='cpu'."
            ) from exc

        return KernelBackend(
            name="cupy",
            xp=cp,
            fft=cupy_fft,
            fftfreq=cupy_fftfreq,
            asnumpy=cp.asnumpy,
        )
    raise ValueError("backend must be 'cpu', 'numpy', 'cupy', or 'gpu'")


def to_numpy(array, backend="cpu"):
    """Move a backend array to a NumPy array."""

    return resolve_backend(backend).asnumpy(array)


def as_blocks(n_frames, num_blocks):
    """Split a trajectory length into equally sized contiguous blocks."""

    if num_blocks < 1:
        raise ValueError("num_blocks must be at least 1")
    block_size = n_frames // num_blocks
    if block_size < 2:
        raise ValueError("each block must contain at least two frames")
    return [(i * block_size, (i + 1) * block_size) for i in range(num_blocks)]


def positive_frequency_slice(num_frames):
    """Return the positive-frequency half used by pySED spectra."""

    return slice(0, num_frames // 2)


def frequency_grid(block_size, dt, positive_only=True, backend="cpu"):
    """Return THz frequencies and the slice used for spectral outputs."""

    kernel = resolve_backend(backend)
    freq = kernel.fftfreq(block_size, dt) / 1e12
    freq_slice = positive_frequency_slice(block_size) if positive_only else slice(None)
    return FrequencyGrid(frequencies_thz=to_numpy(freq[freq_slice], kernel), freq_slice=freq_slice)


def standard_error(samples, backend="cpu"):
    """Return the standard error over the first axis."""

    kernel = resolve_backend(backend)
    xp = kernel.xp
    samples = xp.asarray(samples, dtype=float)
    if samples.shape[0] < 2:
        return np.zeros(samples.shape[1:], dtype=float)
    sem = xp.std(samples, axis=0, ddof=1) / xp.sqrt(samples.shape[0])
    return to_numpy(sem, kernel)


def phase_factors(positions_block, q_vec, backend="cpu"):
    """Return ``exp(i Q dot r_i(t))`` for one trajectory block and Q vector."""

    kernel = resolve_backend(backend)
    xp = kernel.xp
    block = xp.asarray(positions_block, dtype=float)
    q_vec = xp.asarray(q_vec, dtype=float)
    if block.ndim != 3 or block.shape[2] != 3:
        raise ValueError("positions_block must have shape (num_frames, num_atoms, 3)")
    if q_vec.shape != (3,):
        raise ValueError("q_vec must have shape (3,)")
    return xp.exp(1j * xp.tensordot(block, q_vec, axes=([2], [0])))


def density_amplitude(positions_block, q_vec, weights=None, backend="cpu"):
    """Return weighted density amplitude ``rho(Q,t)`` for one block."""

    phase = phase_factors(positions_block, q_vec, backend=backend)
    return density_amplitude_from_phase(phase, weights, backend=backend)


def density_amplitude_from_phase(phase, weights=None, backend="cpu"):
    """Return weighted density amplitude from precomputed phase factors."""

    kernel = resolve_backend(backend)
    xp = kernel.xp
    phase = xp.asarray(phase)
    if phase.ndim != 2:
        raise ValueError("phase must have shape (num_frames, num_atoms)")
    if weights is None:
        weights = xp.ones(phase.shape[1], dtype=float)
    weights = xp.asarray(weights)
    if weights.shape != (phase.shape[1],):
        raise ValueError("weights length must match number of atoms")
    return phase @ weights


def fft_power(signal, axis=0, norm="forward", backend="cpu"):
    """Return ``abs(FFT(signal))**2`` with the configured FFT normalization."""

    kernel = resolve_backend(backend)
    xp = kernel.xp
    signal = xp.asarray(signal)
    return xp.abs(kernel.fft(signal, axis=axis, norm=norm)) ** 2


def circular_autocorrelation(signal, backend="cpu"):
    """Return circular autocorrelation used for finite-block reference DSF."""

    kernel = resolve_backend(backend)
    xp = kernel.xp
    signal = xp.asarray(signal)
    if signal.ndim != 1:
        raise ValueError("signal must be one-dimensional")
    corr = xp.empty(signal.size, dtype=complex)
    signal_conj = xp.conjugate(signal)
    for lag in range(signal.size):
        corr[lag] = xp.mean(xp.roll(signal, -lag) * signal_conj)
    return corr

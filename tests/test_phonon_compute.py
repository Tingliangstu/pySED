"""Numerical and scheduling tests for the SED computation kernel."""

import importlib.util
import unittest
from contextlib import redirect_stdout
from functools import partial
from io import StringIO
from types import SimpleNamespace

import numpy as np
from scipy.fftpack import fft

from pySED.sed.Phonon import _element_symbol_from_mass, spectral_energy_density
from pySED.io.schema import TrajectoryBatch, POSITIONS, VELOCITIES
from pySED.sed.sed_kernel import (
    JOULE_PER_EV,
    JOULE_SECOND_TO_EV_PER_THZ,
    calculate_q_batch,
    finalize_spectrum,
    SEDKernelConfig,
)
from pySED.backends.gpu import (
    GPUBackend,
    _GPUEventProfiler,
    _load_cupy,
)
from pySED.backends.common import (
    qpoint_batches,
)


def _reference_qpoints(vels, cell_vecs, config):
    """Reproduce the previous one-Q-at-a-time implementation."""

    num_frames = vels.shape[0]
    num_qpoints = config.qpoints.shape[0]
    if config.output_partial:
        result = np.zeros(
            (num_frames, num_qpoints, config.num_types, 3)
        )
    else:
        result = np.zeros((num_frames, num_qpoints))

    for q_index, qpoint in enumerate(config.qpoints):
        phase = np.exp(1.0j * np.dot(cell_vecs, qpoint))
        for basis_index, atom_ids in enumerate(config.basis_ids):
            projected_vels = np.tensordot(
                vels[:, atom_ids, :],
                phase,
                axes=([1], [0]),
            )
            spectrum = (
                np.abs(fft(projected_vels, axis=0)) ** 2
                * config.masses[basis_index]
                * config.amu_2_kg
            )
            if config.output_partial:
                type_index = config.basis_to_type[basis_index]
                result[:, q_index, type_index, :] += spectrum
            else:
                result[:, q_index] += spectrum.sum(axis=1)
    return result


class SEDKernelTests(unittest.TestCase):
    """Check that batching and backend changes preserve SED values."""

    @classmethod
    def setUpClass(cls):
        """Create deterministic velocity, topology, and Q-point fixtures."""

        generator = np.random.default_rng(20260726)
        cls.vels = generator.normal(size=(8, 6, 3))
        cls.cell_vecs = generator.normal(size=(3, 3))
        cls.qpoints = generator.normal(size=(7, 3))
        cls.basis_ids = (
            np.array([0, 2, 4]),
            np.array([1, 3, 5]),
        )

    def _config(self, output_partial):
        """Build the small immutable config used by all execution paths."""

        return SEDKernelConfig(
            basis_ids=self.basis_ids,
            masses=np.array([28.085, 15.999]),
            basis_to_type=(1, 0),
            qpoints=self.qpoints,
            num_types=2,
            output_partial=output_partial,
            amu_2_kg=1.66054e-27,
        )

    def test_batched_kernel_matches_previous_total_formula(self):
        """The batched total SED must match the old Q-by-Q loop."""

        config = self._config(output_partial=False)
        actual = calculate_q_batch(
            self.vels,
            self.cell_vecs,
            config.qpoints,
            config.basis_ids,
            config,
        )
        expected = _reference_qpoints(
            self.vels,
            self.cell_vecs,
            config,
        )
        np.testing.assert_allclose(actual, expected, rtol=1e-13)

    def test_batched_kernel_preserves_types_and_xyz_directions(self):
        """Partial output must retain independent atom-type and xyz axes."""

        config = self._config(output_partial=True)
        actual = calculate_q_batch(
            self.vels,
            self.cell_vecs,
            config.qpoints,
            config.basis_ids,
            config,
        )
        expected = _reference_qpoints(
            self.vels,
            self.cell_vecs,
            config,
        )
        self.assertEqual(actual.shape, (8, 7, 2, 3))
        np.testing.assert_allclose(actual, expected, rtol=1e-13)

    def test_qpoint_batches_cover_each_qpoint_once(self):
        """Small worker batches must cover 20-40 Q points without overlap."""

        for num_qpoints in (20, 40):
            for num_workers in (4, 36):
                batches = qpoint_batches(num_qpoints, num_workers)
                flattened = [
                    q_index
                    for start, stop in batches
                    for q_index in range(start, stop)
                ]
                self.assertEqual(flattened, list(range(num_qpoints)))
                self.assertTrue(
                    all(
                        2 <= stop - start <= 5
                        for start, stop in batches
                    )
                )

    def test_gpu_q_batches_make_one_result_transfer_per_block(self):
        """Split or full GPU batches must return one complete block result."""

        class NumPyGPUAdapter:
            """Provide the small CuPy interface used by the GPU method."""

            float64 = np.float64

            def __init__(self):
                self.fft = SimpleNamespace(fft=fft)
                self.transfer_count = 0
                self.cuda = SimpleNamespace(
                    runtime=SimpleNamespace(memGetInfo=lambda: (10**12, 10**12)),
                )

            def __getattr__(self, name):
                return getattr(np, name)

            def asnumpy(self, values):
                self.transfer_count += 1
                return np.asarray(values)

        config = self._config(output_partial=True)
        expected = calculate_q_batch(
            self.vels,
            self.cell_vecs,
            config.qpoints,
            config.basis_ids,
            config,
        )

        for batch_size in (3, len(config.qpoints)):
            adapter = NumPyGPUAdapter()
            calculator = GPUBackend(config, config.qpoints)
            calculator._cupy = adapter
            calculator._gpu_basis_ids = config.basis_ids
            calculator._gpu_qpoints = config.qpoints
            calculator._gpu_q_batch_size = batch_size
            calculator._gpu_event_profiler = SimpleNamespace(
                start=lambda: None,
                stop=lambda *args: None,
                measure=lambda stage, operation, *args, **kwargs: operation(*args, **kwargs),
            )
            result = np.zeros_like(expected)
            with redirect_stdout(StringIO()):
                calculator.compute(self.vels, self.cell_vecs, result)

            self.assertEqual(adapter.transfer_count, 1)
            np.testing.assert_allclose(
                result,
                expected,
                rtol=1e-13,
            )

    def test_gpu_memory_prefers_all_qpoints_when_they_fit(self):
        """The GPU should automatically select one full Q-point batch."""

        config = self._config(output_partial=True)
        runtime = SimpleNamespace(
            memGetInfo=lambda: (10**12, 10**12),
        )
        calculator = GPUBackend(config, config.qpoints)
        calculator._cupy = SimpleNamespace(
            cuda=SimpleNamespace(runtime=runtime),
        )
        result = np.zeros((8, 7, 2, 3))

        with redirect_stdout(StringIO()):
            calculator._check_gpu_memory(
                self.vels,
                self.cell_vecs,
                result,
            )

        self.assertEqual(calculator._gpu_q_batch_size, 7)

    def test_element_symbol_is_inferred_only_for_matching_mass(self):
        """Mass labels should identify standard elements without guessing."""

        self.assertEqual(_element_symbol_from_mass(15.9994), "O")
        self.assertEqual(_element_symbol_from_mass(123.456), "unknown")

    def test_online_blocks_and_one_sided_fold_preserve_power(self):
        """Block averaging and frequency folding must preserve total power."""

        rng = np.random.default_rng(2468)
        stored_blocks = rng.random((3, 8, 7))

        block_average = stored_blocks.sum(axis=0) / stored_blocks.shape[0]
        expected = block_average.copy()
        n_half = expected.shape[0] // 2
        negative = block_average[:n_half:-1, ...]
        paired_length = negative.shape[0]
        expected[1 : 1 + paired_length, ...] = (
            block_average[1 : 1 + paired_length, ...] + negative
        )
        expected = expected[:n_half, ...]
        expected *= JOULE_SECOND_TO_EV_PER_THZ
        retained_two_sided_power = (
            block_average[0, ...]
            + block_average[1 : 1 + paired_length, ...].sum(axis=0)
            + negative.sum(axis=0)
        ) * JOULE_SECOND_TO_EV_PER_THZ

        calculator = object.__new__(spectral_energy_density)
        calculator._sed_sum = np.zeros_like(stored_blocks[0])
        for block in stored_blocks:
            calculator._sed_sum += block
        calculator.num_blocks = stored_blocks.shape[0]
        calculator.freq_fft = np.arange(stored_blocks.shape[1], dtype=float)

        calculator._average_blocks_and_frequencies()

        np.testing.assert_allclose(calculator.sed_avg, expected)
        np.testing.assert_allclose(
            calculator.sed_avg.sum(axis=0),
            retained_two_sided_power,
        )
        self.assertFalse(hasattr(calculator, "_sed_sum"))

    def test_one_sided_fold_handles_odd_and_even_frame_counts(self):
        """Retain every positive bin, excluding only an even-length Nyquist."""
        cases = ((1, [1]), (2, [1]), (3, [1, 5]), (4, [1, 6]),
                 (5, [1, 7, 7]), (6, [1, 8, 8]))
        for num_frames, expected_power in cases:
            for shape in ((2,), (2, 2, 3)):
                with self.subTest(frames=num_frames, shape=shape):
                    scales = np.arange(1, np.prod(shape) + 1).reshape(shape)
                    power = np.arange(1, num_frames + 1, dtype=float)
                    spectrum = power.reshape((-1,) + (1,) * len(shape)) * scales
                    frequencies = np.fft.fftfreq(num_frames, d=0.001)
                    actual, actual_frequencies = finalize_spectrum(
                        2 * spectrum, frequencies, num_blocks=2,
                    )
                    expected = np.asarray(expected_power).reshape(
                        (-1,) + (1,) * len(shape)
                    ) * scales * JOULE_SECOND_TO_EV_PER_THZ
                    np.testing.assert_allclose(actual, expected)
                    np.testing.assert_array_equal(
                        actual_frequencies, frequencies[frequencies >= 0],
                    )

    def test_five_frame_signal_retains_highest_positive_frequency_peak(self):
        """A 400 THz signal sampled at 1 fs must survive the final fold."""
        velocities = np.zeros((5, 1, 3))
        velocities[:, 0, 0] = np.cos(2 * np.pi * 2 * np.arange(5) / 5)
        config = SEDKernelConfig(
            basis_ids=(np.array([0]),), masses=np.array([1.0]),
            basis_to_type=(0,), qpoints=np.zeros((1, 3)), num_types=1,
            output_partial=False, amu_2_kg=1.0,
        )
        spectrum = calculate_q_batch(
            velocities, np.zeros((1, 3)), config.qpoints,
            config.basis_ids, config,
        )
        actual, frequencies = finalize_spectrum(
            spectrum, np.fft.fftfreq(5, d=0.001), num_blocks=1,
        )
        np.testing.assert_array_equal(frequencies, [0.0, 200.0, 400.0])
        np.testing.assert_allclose(
            actual[:, 0] / JOULE_SECOND_TO_EV_PER_THZ,
            [0.0, 0.0, 12.5], atol=1e-12,
        )

    def test_output_unit_conversion_preserves_integrated_energy(self):
        """The eV/THz spectrum must integrate to the same physical energy."""

        frequency_thz = np.linspace(0.0, 10.0, 101)
        sed_js = 2.0e-33 * (1.0 + frequency_thz)
        sed_ev_per_thz = (
            sed_js * JOULE_SECOND_TO_EV_PER_THZ
        )

        energy_from_internal_ev = (
            2
            * np.pi
            * 1.0e12
            * np.trapezoid(sed_js, frequency_thz)
            / JOULE_PER_EV
        )
        energy_from_output_ev = np.trapezoid(
            sed_ev_per_thz,
            frequency_thz,
        )

        self.assertAlmostEqual(
            energy_from_output_ev,
            energy_from_internal_ev,
            places=14,
        )

    def test_driver_legacy_scale_matches_2x_normalization(self):
        """Check both scales against the 2.x FFT formula, including DC and odd N."""
        for num_frames in (5, 8):
            velocities = self.vels[:num_frames]
            positions = np.broadcast_to(
                np.repeat(self.cell_vecs, 2, axis=0), velocities.shape,
            )
            blocks = [
                TrajectoryBatch(start=i * num_frames, data={
                    POSITIONS: positions, VELOCITIES: velocities * amplitude,
                })
                for i, amplitude in enumerate((1.0, 2.0))
            ]
            for partial_output in (False, True):
                config = self._config(partial_output)
                lattice = SimpleNamespace(
                    unitcell_index=np.repeat([1, 2, 3], 2),
                    basis_index=np.tile([1, 2], 3),
                    masses=np.tile(config.masses, 3),
                    cell_ref_ids=np.array([0, 2, 4]), num_qpoints=7,
                    qpoints=config.qpoints, reduced_qpoints=config.qpoints,
                )
                for spacing_fs in (1.0, 2.0):
                    dt = spacing_fs * 1e-15
                    raw = _reference_qpoints(velocities, self.cell_vecs, config)
                    # The two block powers are 1 and 4 times the reference.
                    old_two_sided = raw * 2.5 / (4 * np.pi * num_frames * dt * 3)
                    frequencies = np.fft.fftfreq(num_frames, dt) / 1e12
                    indices = np.flatnonzero(frequencies >= 0)
                    old_folded = old_two_sided[indices].copy()
                    for index in indices[1:]:
                        old_folded[index] = (
                            old_two_sided[index] + old_two_sided[-index]
                        ) / 2.0
                    for legacy in (False, True):
                        with self.subTest(frames=num_frames, partial=partial_output,
                                          dt=dt, legacy=legacy):
                            params = SimpleNamespace(
                                total_num_steps=2 * num_frames, output_data_stride=1,
                                num_blocks=2, num_atoms=6, time_step=spacing_fs,
                                backend="numpy", max_cores=1,
                                output_partial=partial_output, trajectory_prefetch=False,
                                legacy_unit=legacy,
                            )
                            with redirect_stdout(StringIO()):
                                calculator = spectral_energy_density(
                                    params, SimpleNamespace(n_atoms=6),
                                )
                                calculator.compute_sed(
                                    params, lattice, blocks[0], iter(blocks[1:]),
                                )
                            expected = old_folded.copy()
                            if not legacy:
                                expected[1:] *= 2.0
                                expected *= dt**2 * JOULE_SECOND_TO_EV_PER_THZ
                            np.testing.assert_allclose(
                                calculator.sed_avg, expected, rtol=1e-12, atol=0,
                            )
                            np.testing.assert_array_equal(
                                calculator.freq_fft, frequencies[indices],
                            )

    def test_block_scaling_includes_discrete_time_step_squared(self):
        """The FFT sum must include dt twice after taking its magnitude squared."""

        calculator = object.__new__(spectral_energy_density)
        calculator.dt = 2.0
        calculator.t_o = 8.0
        calculator.num_unit_cells = 3
        calculator.num_blocks = 1
        calculator._sed_sum = np.zeros((8, 7))
        calculator.trajectory_prefetch = False
        backend = SimpleNamespace(
            prepare=lambda *args, **kwargs: None,
            compute=lambda vels, cells, result: result.fill(1.0),
        )
        calculator._loop_over_blocks(
            backend=backend,
            trajectory_blocks=[(self.vels, self.cell_vecs)],
        )

        expected = calculator.dt**2 / (
            4
            * np.pi
            * calculator.t_o
            * calculator.num_unit_cells
        )
        np.testing.assert_allclose(calculator._sed_sum, expected)

    @unittest.skipUnless(
        importlib.util.find_spec("cupy") is not None,
        "CuPy is not installed",
    )
    def test_cupy_kernel_matches_numpy_when_cuda_is_available(self):
        """The optional GPU backend must use the same batched SED formula."""

        try:
            cp = _load_cupy()
        except RuntimeError as error:
            self.skipTest(str(error))

        config = self._config(output_partial=True)
        expected = calculate_q_batch(
            self.vels,
            self.cell_vecs,
            config.qpoints,
            config.basis_ids,
            config,
        )
        profiler = _GPUEventProfiler(cp)
        actual = calculate_q_batch(
            cp.asarray(self.vels),
            cp.asarray(self.cell_vecs),
            cp.asarray(config.qpoints),
            tuple(cp.asarray(ids) for ids in config.basis_ids),
            config,
            xp=cp,
            fft_function=partial(profiler.measure, "fft", cp.fft.fft),
            project_function=partial(profiler.measure, "projection", cp.tensordot),
        )
        event_times = profiler.elapsed_seconds()
        np.testing.assert_allclose(
            cp.asnumpy(actual),
            expected,
            rtol=1e-11,
        )
        self.assertGreater(event_times["projection"], 0.0)
        self.assertGreater(event_times["fft"], 0.0)


if __name__ == "__main__":
    unittest.main()

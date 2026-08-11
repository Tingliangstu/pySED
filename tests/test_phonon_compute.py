"""Numerical and scheduling tests for the SED computation kernel."""

import importlib.util
import unittest
from concurrent.futures import ProcessPoolExecutor
from contextlib import redirect_stdout
from io import StringIO
from types import SimpleNamespace

import numpy as np
from scipy.fftpack import fft

from pySED.sed.Phonon import (
    JOULE_PER_EV,
    JOULE_SECOND_TO_EV_PER_THZ,
    _calculate_q_batch,
    _compute_shared_q_batch,
    _element_symbol_from_mass,
    _estimate_kernel_working_bytes,
    _GPUEventProfiler,
    _initialize_cpu_worker,
    _load_cupy,
    _qpoint_batches,
    _SEDKernelConfig,
    _share_numpy_array,
    spectral_energy_density,
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

        return _SEDKernelConfig(
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
        actual = _calculate_q_batch(
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
        actual = _calculate_q_batch(
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
                batches = _qpoint_batches(num_qpoints, num_workers)
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

    def test_qpoint_progress_keeps_original_output_style(self):
        """Internal batching must not change the visible per-Q progress text."""

        calculator = object.__new__(spectral_energy_density)
        calculator.reduced_qpoints = np.array(
            [[0.0, 0.25, 0.5], [0.5, 0.25, 0.0]]
        )
        output = StringIO()
        with redirect_stdout(output):
            calculator._print_qpoints(0, 2, 2)

        self.assertEqual(
            output.getvalue(),
            "\tNow calculating q-point 1/2:\t"
            "q = (0.0000, 0.2500, 0.5000)\n"
            "\tNow calculating q-point 2/2:\t"
            "q = (0.5000, 0.2500, 0.0000)\n",
        )

    def test_gpu_timing_places_other_work_below_the_six_row_sum(self):
        """The complete timing reconciliation should read as one table."""

        calculator = object.__new__(spectral_energy_density)
        calculator._gpu_timings = {"trajectory": 2.0}
        calculator._gpu_event_profiler = SimpleNamespace(
            elapsed_seconds=lambda: {
                "compute": 3.0,
                "projection": 1.0,
                "fft": 0.5,
                "upload": 0.4,
                "download": 0.1,
            }
        )
        calculator._gpu_total_bytes = 0

        output = StringIO()
        with redirect_stdout(output):
            calculator._print_gpu_timing_profile(total_elapsed=6.5)

        lines = output.getvalue().splitlines()
        sum_index = next(
            index
            for index, line in enumerate(lines)
            if "Sum of the six rows above" in line
        )
        self.assertIn("5.50 s", lines[sum_index])
        self.assertIn("Other SED work", lines[sum_index + 1])
        self.assertIn("1.00 s", lines[sum_index + 1])

    def test_gpu_q_batches_make_one_result_transfer_per_block(self):
        """Split or full GPU batches must return one complete block result."""

        class NumPyGPUAdapter:
            """Provide the small CuPy interface used by the GPU method."""

            float64 = np.float64

            def __init__(self):
                self.fft = SimpleNamespace(fft=fft)
                self.transfer_count = 0

            def __getattr__(self, name):
                return getattr(np, name)

            def asnumpy(self, values):
                self.transfer_count += 1
                return np.asarray(values)

        config = self._config(output_partial=True)
        expected = _calculate_q_batch(
            self.vels,
            self.cell_vecs,
            config.qpoints,
            config.basis_ids,
            config,
        )

        for batch_size in (3, len(config.qpoints)):
            adapter = NumPyGPUAdapter()
            calculator = object.__new__(spectral_energy_density)
            calculator._cupy = adapter
            calculator._gpu_basis_ids = config.basis_ids
            calculator._gpu_qpoints = config.qpoints
            calculator._gpu_q_batch_size = batch_size
            calculator.qdot = np.zeros_like(expected)
            calculator._print_qpoints = lambda *args: None

            calculator._compute_qpoints_gpu(
                self.vels,
                self.cell_vecs,
                config,
            )

            self.assertEqual(adapter.transfer_count, 1)
            np.testing.assert_allclose(
                calculator.qdot,
                expected,
                rtol=1e-13,
            )

    def test_gpu_memory_prefers_all_qpoints_when_they_fit(self):
        """The GPU should automatically select one full Q-point batch."""

        config = self._config(output_partial=True)
        runtime = SimpleNamespace(
            memGetInfo=lambda: (10**12, 10**12),
        )
        calculator = object.__new__(spectral_energy_density)
        calculator._cupy = SimpleNamespace(
            cuda=SimpleNamespace(runtime=runtime),
        )
        calculator.qdot = np.zeros((8, 7, 2, 3))

        with redirect_stdout(StringIO()):
            calculator._check_gpu_memory(
                self.vels,
                self.cell_vecs,
                config,
            )

        self.assertEqual(calculator._gpu_q_batch_size, 7)

    def test_element_symbol_is_inferred_only_for_matching_mass(self):
        """Mass labels should identify standard elements without guessing."""

        self.assertEqual(_element_symbol_from_mass(15.9994), "O")
        self.assertEqual(_element_symbol_from_mass(123.456), "unknown")

    def test_cpu_memory_estimate_uses_real_array_shapes(self):
        """The CPU estimate should use block arrays and keep cute output."""

        config = self._config(output_partial=False)
        working_bytes = _estimate_kernel_working_bytes(
            self.vels,
            config,
            max_q_batch=4,
        )
        self.assertGreater(working_bytes, 0)

        calculator = object.__new__(spectral_energy_density)
        calculator.backend = "numpy"
        calculator._sed_sum = np.zeros((8, 7))
        calculator.qdot = np.zeros((8, 7))
        output = StringIO()
        with redirect_stdout(output):
            calculator._print_cpu_memory_estimate(
                self.vels,
                self.cell_vecs,
                config,
                executor=None,
                num_workers=1,
            )
        self.assertIn(
            "GB (trajectory block reading + SED arrays;",
            output.getvalue(),
        )

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

    def test_block_scaling_includes_discrete_time_step_squared(self):
        """The FFT sum must include dt twice after taking its magnitude squared."""

        config = self._config(output_partial=False)
        calculator = object.__new__(spectral_energy_density)
        calculator.dt = 2.0
        calculator.t_o = 8.0
        calculator.num_unit_cells = 3
        calculator.num_blocks = 1
        calculator._sed_sum = np.zeros((8, 7))
        calculator._allocate_qdot = lambda num_qpoints: setattr(
            calculator,
            "qdot",
            np.ones((8, num_qpoints)),
        )
        calculator._print_cpu_memory_estimate = lambda *args: None
        calculator._loop_over_qpoints = lambda *args: None

        lattice = SimpleNamespace(num_qpoints=7)
        calculator._loop_over_blocks(
            lattice_info=lattice,
            kernel_config=config,
            executor=None,
            num_workers=1,
            trajectory_blocks=[(self.vels, self.cell_vecs)],
        )

        expected = calculator.dt**2 / (
            4
            * np.pi
            * calculator.t_o
            * calculator.num_unit_cells
        )
        np.testing.assert_allclose(calculator._sed_sum, expected)

    def test_persistent_workers_read_shared_velocity_blocks(self):
        """Two blocks can reuse one pool without serializing full velocities."""

        config = self._config(output_partial=False)
        executor = ProcessPoolExecutor(
            max_workers=2,
            initializer=_initialize_cpu_worker,
            initargs=(config,),
        )
        try:
            for scale in (1.0, 2.0):
                block_vels = self.vels * scale
                expected = _calculate_q_batch(
                    block_vels,
                    self.cell_vecs,
                    config.qpoints,
                    config.basis_ids,
                    config,
                )
                actual = np.zeros_like(expected)

                with _share_numpy_array(block_vels) as shared:
                    name, shape, dtype_string = shared
                    futures = [
                        executor.submit(
                            _compute_shared_q_batch,
                            name,
                            shape,
                            dtype_string,
                            self.cell_vecs,
                            start,
                            stop,
                        )
                        for start, stop in _qpoint_batches(7, 2)
                    ]
                    for future in futures:
                        start, stop, values = future.result()
                        actual[:, start:stop] = values

                np.testing.assert_allclose(actual, expected, rtol=1e-13)
        finally:
            executor.shutdown()

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
        expected = _calculate_q_batch(
            self.vels,
            self.cell_vecs,
            config.qpoints,
            config.basis_ids,
            config,
        )
        profiler = _GPUEventProfiler(cp)
        actual = _calculate_q_batch(
            cp.asarray(self.vels),
            cp.asarray(self.cell_vecs),
            cp.asarray(config.qpoints),
            tuple(cp.asarray(ids) for ids in config.basis_ids),
            config,
            xp=cp,
            fft_function=cp.fft.fft,
            gpu_profiler=profiler,
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

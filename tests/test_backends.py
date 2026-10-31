"""Resource ownership and import boundaries of the split SED implementation."""

import subprocess
import sys
import unittest
from contextlib import redirect_stdout
from dataclasses import replace
from io import StringIO
from multiprocessing import shared_memory
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np

from pySED.backends.cpu import CPUBackend, _share_numpy_array
from pySED.backends.gpu import GPUBackend
from pySED.sed.Phonon import spectral_energy_density
from pySED.io.schema import POSITIONS, VELOCITIES, TrajectoryBatch
from pySED.sed.sed_kernel import SEDKernelConfig, calculate_q_batch


class BackendTests(unittest.TestCase):
    def setUp(self):
        rng = np.random.default_rng(20261030)
        self.vels = rng.normal(size=(8, 6, 3))
        self.positions = rng.normal(size=(8, 6, 3))
        self.cell_vecs = self.positions[:, [0, 2, 4], :].mean(axis=0)
        self.config = SEDKernelConfig(
            basis_ids=(np.array([0, 2, 4]), np.array([1, 3, 5])),
            masses=np.array([28.085, 15.999]),
            basis_to_type=(1, 0),
            qpoints=rng.normal(size=(7, 3)),
            num_types=2,
            output_partial=False,
            amu_2_kg=1.66054e-27,
        )

    def test_serial_and_persistent_cpu_backends_match_for_two_blocks(self):
        for partial_output in (False, True):
            config = replace(self.config, output_partial=partial_output)
            for cores in (1, 2):
                with self.subTest(partial=partial_output, cores=cores):
                    with redirect_stdout(StringIO()):
                        with CPUBackend(config, config.qpoints, cores) as backend:
                            pool = backend._executor
                            for scale in (1.0, 2.0):
                                vels = self.vels * scale
                                expected = calculate_q_batch(
                                    vels, self.cell_vecs, config.qpoints,
                                    config.basis_ids, config,
                                )
                                result = np.empty_like(expected)
                                backend.compute(vels, self.cell_vecs, result)
                                np.testing.assert_allclose(result, expected, rtol=1e-13)
                                self.assertIs(backend._executor, pool)
                    self.assertIsNone(backend._executor)

    def test_worker_error_propagates_and_closes_the_pool(self):
        config = replace(
            self.config,
            basis_ids=(np.array([0, 2, 999]), self.config.basis_ids[1]),
        )
        backend = CPUBackend(config, config.qpoints, 2)
        with redirect_stdout(StringIO()):
            with self.assertRaises(IndexError):
                with backend:
                    backend.compute(self.vels, self.cell_vecs, np.zeros((8, 7)))
        self.assertIsNone(backend._executor)

    def test_shared_buffer_is_removed_after_an_exception(self):
        with self.assertRaisesRegex(RuntimeError, "deliberate failure"):
            with _share_numpy_array(self.vels) as (name, shape, dtype):
                self.assertEqual(shape, self.vels.shape)
                raise RuntimeError("deliberate failure")
        with self.assertRaises(FileNotFoundError):
            shared_memory.SharedMemory(name=name)

    def test_gpu_initialization_failure_releases_owned_buffers(self):
        device_pool, pinned_pool = Mock(), Mock()
        cp = SimpleNamespace(
            asarray=Mock(side_effect=RuntimeError("upload failed")),
            get_default_memory_pool=lambda: device_pool,
            get_default_pinned_memory_pool=lambda: pinned_pool,
        )
        backend = GPUBackend(self.config, self.config.qpoints)
        with patch("pySED.backends.gpu._load_cupy", return_value=cp):
            with self.assertRaisesRegex(RuntimeError, "upload failed"):
                with backend:
                    self.fail("initialization must fail before entering")
        device_pool.free_all_blocks.assert_called_once()
        pinned_pool.free_all_blocks.assert_called_once()
        self.assertIsNone(backend._gpu_basis_ids)
        self.assertIsNone(backend._gpu_qpoints)

    def test_driver_serial_and_parallel_outputs_match(self):
        lattice = SimpleNamespace(
            unitcell_index=np.repeat([1, 2, 3], 2),
            basis_index=np.tile([1, 2], 3),
            masses=np.tile(self.config.masses, 3),
            cell_ref_ids=np.array([0, 2, 4]),
            num_qpoints=7,
            qpoints=self.config.qpoints,
            reduced_qpoints=self.config.qpoints,
        )
        blocks = [
            TrajectoryBatch(
                start=index * 8,
                data={POSITIONS: self.positions, VELOCITIES: self.vels * scale},
            )
            for index, scale in enumerate((1.0, 2.0))
        ]
        for partial_output in (False, True):
            results = []
            for cores in (1, 2):
                params = SimpleNamespace(
                    total_num_steps=16, output_data_stride=1, num_blocks=2,
                    num_atoms=6, time_step=0.5, backend="numpy", max_cores=cores,
                    output_partial=partial_output, trajectory_prefetch=False,
                )
                with redirect_stdout(StringIO()):
                    calculator = spectral_energy_density(params, SimpleNamespace(n_atoms=6))
                    calculator.compute_sed(params, lattice, blocks[0], iter(blocks[1:]))
                results.append(calculator.sed_avg)
            np.testing.assert_allclose(results[0], results[1], rtol=1e-13)

    def test_cpu_import_does_not_load_cupy(self):
        command = (
            "import sys; from pySED.sed import sed_kernel; "
            "assert 'pySED.backends' not in sys.modules; "
            "assert 'pySED.io' not in sys.modules; "
            "from pySED.postprocess import Lorentz; "
            "assert 'matplotlib' not in sys.modules; "
            "assert 'pySED.postprocess.FileIO' not in sys.modules; "
            "from pySED.sed import Phonon; "
            "from pySED.postprocess import FileIO, Lorentz, Plot_SED; "
            "assert 'cupy' not in sys.modules; "
            "assert 'pySED.backends.gpu' not in sys.modules"
        )
        subprocess.run(
            [sys.executable, "-c", command],
            cwd=Path(__file__).resolve().parents[1], check=True,
            capture_output=True, text=True,
        )

if __name__ == "__main__":
    unittest.main()

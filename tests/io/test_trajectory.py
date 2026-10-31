"""Tests for direct/cached trajectory blocks and one-block prefetch."""

import importlib.util
import threading
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import numpy as np

from pySED.io.prepare import prepare_trajectory
from pySED.io.schema import CELLS, POSITIONS, VELOCITIES
from pySED.io.text import GPUMD_XYZ, LAMMPS_DUMP
from pySED.io.trajectory import (
    DirectTextBlockSource,
    NetCDFBlockSource,
    prefetch_one,
)
from pySED.sed.Phonon import spectral_energy_density
from tests.io.test_text import GPUMD_TEXT, LAMMPS_TEXT

NETCDF4_AVAILABLE = importlib.util.find_spec("netCDF4") is not None


def _gpumd_frames(count=8):
    frames = []
    for index in range(count):
        shift = 0.01 * index
        frames.extend(
            [
                "2",
                (
                    'Lattice="10 0 0 0 10 0 0 0 10" '
                    "Properties=species:S:1:pos:R:3:vel:R:3"
                ),
                f"Sr {shift} 0 0 {0.001 + shift} 0.002 0.003",
                f"Ti {1 + shift} 1 1 0.004 {0.005 + shift} 0.006",
            ]
        )
    return "\n".join(frames) + "\n"


class DirectTrajectoryTests(unittest.TestCase):
    def test_gpumd_and_lammps_use_the_same_block_interface(self):
        cases = (
            ("dump.xyz", GPUMD_TEXT, GPUMD_XYZ, "metal"),
            ("dump.lammpstrj", LAMMPS_TEXT, LAMMPS_DUMP, "metal"),
        )
        with TemporaryDirectory() as directory:
            for name, content, source_format, lammps_unit in cases:
                with self.subTest(source_format=source_format):
                    path = Path(directory) / name
                    path.write_text(content, encoding="utf-8")
                    source = DirectTextBlockSource(
                        path=path,
                        source_format=source_format,
                        required_frames=2,
                        parser_batch_size=1,
                        lammps_unit=lammps_unit,
                    )
                    blocks = list(source.iter_blocks(1, 2))

                    self.assertEqual([block.start for block in blocks], [0, 1])
                    self.assertEqual(source.info.n_atoms, 2)
                    for block in blocks:
                        self.assertEqual(
                            set(block.data),
                            {POSITIONS, VELOCITIES, CELLS},
                        )
                        self.assertEqual(block.data[POSITIONS].shape, (1, 2, 3))

                    with self.assertRaisesRegex(RuntimeError, "consumed once"):
                        next(source.iter_blocks(1, 2))

    def test_prefetch_runs_the_next_item_in_one_background_thread(self):
        thread_names = []

        def items():
            for value in range(3):
                thread_names.append(threading.current_thread().name)
                yield value

        with prefetch_one(items(), enabled=True) as prefetched:
            self.assertEqual(list(prefetched), [0, 1, 2])

        self.assertTrue(thread_names)
        self.assertTrue(
            all(name.startswith("pySED-trajectory") for name in thread_names)
        )


@unittest.skipUnless(NETCDF4_AVAILABLE, "netCDF4 is not installed")
class CachedTrajectoryParityTests(unittest.TestCase):
    def _calculate_sed(self, source, prefetch):
        raw_blocks = source.iter_blocks(block_size=4, num_blocks=2)
        params = SimpleNamespace(
            total_num_steps=8,
            output_data_stride=1,
            num_blocks=2,
            num_atoms=2,
            time_step=1.0,
            device="cpu",
            max_cores=1,
            output_partial=False,
            trajectory_prefetch=prefetch,
        )
        lattice = SimpleNamespace(
            unitcell_index=np.array([1, 1]),
            basis_index=np.array([1, 2]),
            masses=np.array([87.62, 47.867]),
            qpoints=np.array([[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]),
            reduced_qpoints=np.array(
                [[0.0, 0.0, 0.0], [0.1, 0.0, 0.0]]
            ),
            num_qpoints=2,
            cell_ref_ids=np.array([0]),
        )
        with prefetch_one(raw_blocks, enabled=prefetch) as blocks:
            first_block = next(blocks)
            calculator = spectral_energy_density(params, source.info)
            with redirect_stdout(StringIO()):
                calculator.compute_sed(
                    params,
                    lattice,
                    first_block,
                    blocks,
                )
        return calculator.freq_fft, calculator.sed_avg

    def test_direct_and_cached_netcdf_produce_the_same_sed(self):
        with TemporaryDirectory() as directory:
            path = Path(directory) / "dump.xyz"
            path.write_text(_gpumd_frames(), encoding="utf-8")
            cached_path, _ = prepare_trajectory(
                path,
                batch_size=3,
                compression_level=0,
            )

            direct = DirectTextBlockSource(
                path=path,
                source_format=GPUMD_XYZ,
                required_frames=8,
                parser_batch_size=3,
                lammps_unit="metal",
            )
            direct_frequency, direct_sed = self._calculate_sed(
                direct,
                prefetch=False,
            )

            cached = NetCDFBlockSource(cached_path, required_frames=8)
            cached_frequency, cached_sed = self._calculate_sed(
                cached,
                prefetch=True,
            )

            np.testing.assert_array_equal(cached_frequency, direct_frequency)
            np.testing.assert_allclose(
                cached_sed,
                direct_sed,
                rtol=2.0e-6,
                atol=0.0,
            )

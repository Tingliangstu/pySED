"""Tests for the one canonical NetCDF reader and writer API."""

import importlib.util
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from pySED.io.netcdf import NetCDFReader, is_netcdf
from pySED.io.prepare import prepare_trajectory
from pySED.io.schema import ATOM_IDS, POSITIONS, TYPES, VELOCITIES

NETCDF4_AVAILABLE = importlib.util.find_spec("netCDF4") is not None


@unittest.skipUnless(NETCDF4_AVAILABLE, "netCDF4 is not installed")
class NetCDFTests(unittest.TestCase):
    def _create_trajectory(self, path: Path) -> None:
        from netCDF4 import Dataset

        with Dataset(str(path), "w", format="NETCDF3_64BIT_OFFSET") as dataset:
            dataset.program = "GPUMD"
            dataset.createDimension("frame", 2)
            dataset.createDimension("atom", 2)
            dataset.createDimension("spatial", 3)

            coordinates = dataset.createVariable(
                "coordinates",
                "f4",
                ("frame", "atom", "spatial"),
            )
            coordinates.units = "angstrom"
            coordinates[:] = np.arange(12, dtype=np.float32).reshape(2, 2, 3)

            velocities = dataset.createVariable(
                "velocities",
                "f4",
                ("frame", "atom", "spatial"),
            )
            velocities.units = "angstrom/picosecond"
            velocities[:] = np.ones((2, 2, 3), dtype=np.float32)

            atom_ids = dataset.createVariable("atom_id", "i4", ("atom",))
            atom_ids[:] = [10, 20]

            atom_types = dataset.createVariable(
                "type",
                "i4",
                ("frame", "atom"),
            )
            atom_types[:] = [[0, 1], [0, 1]]

            lengths = dataset.createVariable(
                "cell_lengths",
                "f8",
                ("frame", "spatial"),
            )
            lengths.units = "angstrom"
            lengths[:] = [[10, 20, 30], [10, 20, 30]]

            angles = dataset.createVariable(
                "cell_angles",
                "f8",
                ("frame", "spatial"),
            )
            angles.units = "degree"
            angles[:] = [[90, 90, 90], [90, 90, 90]]

    def test_direct_read_normalizes_units_and_supports_slices(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "trajectory.nc"
            self._create_trajectory(source)

            self.assertTrue(is_netcdf(source))
            prepared, source_format = prepare_trajectory(source)
            self.assertEqual(source_format, "gpumd_netcdf")
            with NetCDFReader(prepared) as trajectory:
                self.assertEqual(trajectory.info.n_frames, 2)
                self.assertEqual(trajectory.info.n_atoms, 2)
                self.assertEqual(trajectory.info.source_program, "GPUMD")
                positions = trajectory.read(POSITIONS, slice(1, 2))
                velocities = trajectory.read(VELOCITIES, slice(0, 1))
                atom_types = trajectory.read(TYPES, slice(1, 2))
                atom_ids = trajectory.read(ATOM_IDS, slice(1, 2))
                cells = trajectory.read_cells()

            self.assertEqual(positions.shape, (1, 2, 3))
            np.testing.assert_allclose(velocities, 100.0)
            np.testing.assert_array_equal(atom_types, [[0, 1]])
            np.testing.assert_array_equal(atom_ids, [10, 20])
            np.testing.assert_allclose(
                np.diagonal(cells, axis1=1, axis2=2),
                [[10, 20, 30], [10, 20, 30]],
                atol=1.0e-12,
            )

    def test_batch_reader_has_bounded_shapes(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "trajectory.nc"
            self._create_trajectory(source)

            prepared, _ = prepare_trajectory(source)
            with NetCDFReader(prepared) as trajectory:
                batches = list(
                    trajectory.iter_batches(
                        fields=[POSITIONS, VELOCITIES],
                        batch_size=1,
                    )
                )

            self.assertEqual(len(batches), 2)
            self.assertEqual(batches[0][POSITIONS].shape, (1, 2, 3))


if __name__ == "__main__":
    unittest.main()

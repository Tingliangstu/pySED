"""Tests for fast GPUMD and LAMMPS text-to-NetCDF conversion."""

import importlib.util
import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from pySED.io.netcdf import NetCDFReader
from pySED.io.prepare import default_converted_path, prepare_trajectory
from pySED.io.schema import ATOM_IDS, POSITIONS, TYPES, VELOCITIES
from pySED.io.text import (
    GPUMD_XYZ,
    LAMMPS_DUMP,
    detect_text_format,
    iter_gpumd_xyz,
)

NETCDF4_AVAILABLE = importlib.util.find_spec("netCDF4") is not None


GPUMD_TEXT = """\
2
Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3:vel:R:3
Sr 0 0 0 0.001 0.002 0.003
Ti 1 1 1 0.004 0.005 0.006
2
Lattice="10 0 0 0 10 0 0 0 10" Properties=species:S:1:pos:R:3:vel:R:3
Sr 0.1 0 0 0.002 0.003 0.004
Ti 1.1 1 1 0.005 0.006 0.007
"""


LAMMPS_TEXT = """\
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id type x y z vx vy vz
2 2 2 2 2 4 5 6
1 1 1 1 1 1 2 3
ITEM: TIMESTEP
10
ITEM: NUMBER OF ATOMS
2
ITEM: BOX BOUNDS pp pp pp
0 10
0 20
0 30
ITEM: ATOMS id type x y z vx vy vz
1 1 1.1 1 1 2 3 4
2 2 2.1 2 2 5 6 7
"""


class GPUMDPropertiesTests(unittest.TestCase):
    def test_extxyz_uses_the_declared_properties_order(self) -> None:
        text = """\
2
Lattice="10 0 0 0 10 0 0 0 10" Properties=vel:R:3:species:S:1:pos:R:3
0.001 0.002 0.003 Sr 1 2 3
0.004 0.005 0.006 Ti 4 5 6
"""
        with TemporaryDirectory() as directory:
            source = Path(directory) / "reordered.xyz"
            source.write_text(text, encoding="utf-8")

            batches = list(iter_gpumd_xyz(source, batch_size=1))

        self.assertEqual(len(batches), 1)
        batch = batches[0].data
        np.testing.assert_allclose(batch[POSITIONS], [[[1, 2, 3], [4, 5, 6]]])
        np.testing.assert_allclose(
            batch[VELOCITIES],
            [[[100, 200, 300], [400, 500, 600]]],
        )
        np.testing.assert_array_equal(batch[TYPES], [[0, 1]])


@unittest.skipUnless(NETCDF4_AVAILABLE, "netCDF4 is not installed")
class TextConversionTests(unittest.TestCase):
    def test_float64_conversion_preserves_position_and_velocity_precision(self) -> None:
        position = 1.123456789012345
        velocity = 0.001234567890123456
        xyz_text = GPUMD_TEXT.replace("0.005", f"{velocity:.17g}")
        lammps_text = LAMMPS_TEXT.replace(
            "1.1 1 1 2", f"1.1 1 1 {velocity:.17g}"
        )
        cases = (
            ("dump.xyz", xyz_text, 1, "metal", 100000.0),
            ("dump.data", lammps_text, 0, "metal", 100.0),
            ("dump.data", lammps_text, 0, "real", 100000.0),
        )
        for filename, text, atom, unit, velocity_scale in cases:
            with (
                self.subTest(filename=filename, unit=unit),
                TemporaryDirectory() as directory,
            ):
                source = Path(directory) / filename
                source.write_text(
                    text.replace("1.1 ", f"{position:.17g} "),
                    encoding="utf-8",
                )
                prepared, _ = prepare_trajectory(
                    source, storage_dtype="float64", lammps_unit=unit, batch_size=1
                )
                with NetCDFReader(prepared) as trajectory:
                    positions = trajectory.read(POSITIONS)
                    velocities = trajectory.read(VELOCITIES)
                self.assertEqual(positions.dtype, np.dtype("float64"))
                self.assertEqual(velocities.dtype, np.dtype("float64"))
                self.assertEqual(float(positions[1, atom, 0]), position)
                self.assertEqual(
                    float(velocities[1, atom, 0]), velocity * velocity_scale
                )

    def test_float64_conversion_refreshes_old_or_lower_precision_cache(self) -> None:
        from netCDF4 import Dataset

        value = 1.123456789012345
        cases = (
            ("dump.xyz", GPUMD_TEXT, 1),
            ("dump.data", LAMMPS_TEXT, 0),
        )
        for filename, text, atom in cases:
            for old_reader in (False, True):
                with (
                    self.subTest(filename=filename, old_reader=old_reader),
                    TemporaryDirectory() as directory,
                ):
                    source = Path(directory) / filename
                    source.write_text(
                        text.replace("1.1 ", f"{value:.17g} "), encoding="utf-8"
                    )
                    initial_dtype = "float64" if old_reader else "float32"
                    prepared, _ = prepare_trajectory(
                        source, storage_dtype=initial_dtype
                    )
                    if old_reader:
                        # Version 1 rounded text through float32 even for float64 caches.
                        with Dataset(str(prepared), "a") as dataset:
                            dataset.pysed_schema_version = 1
                            dataset.variables["coordinates"][1, atom, 0] = float(
                                np.float32(value)
                            )
                    with NetCDFReader(prepared) as trajectory:
                        self.assertNotEqual(
                            float(trajectory.read(POSITIONS)[1, atom, 0]), value
                        )

                    refreshed, _ = prepare_trajectory(source, storage_dtype="float64")
                    self.assertEqual(refreshed, prepared)
                    with NetCDFReader(refreshed) as trajectory:
                        self.assertEqual(
                            float(trajectory.read(POSITIONS)[1, atom, 0]), value
                        )

                    first_mtime = refreshed.stat().st_mtime_ns
                    prepare_trajectory(source, storage_dtype="float64")
                    self.assertEqual(refreshed.stat().st_mtime_ns, first_mtime)

    def test_conversion_prints_progress_and_destination(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.xyz"
            source.write_text(GPUMD_TEXT, encoding="utf-8")
            output = io.StringIO()

            with redirect_stdout(output):
                prepared, _ = prepare_trajectory(source, batch_size=1)

            message = output.getvalue()
            self.assertIn("Converting text trajectory to NetCDF", message)
            self.assertIn(str(prepared), message)
            self.assertIn("100%", message)
            self.assertIn("2 frames", message)

    def test_gpumd_xyz_is_streamed_and_reused(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.xyz"
            source.write_text(GPUMD_TEXT, encoding="utf-8")

            self.assertEqual(detect_text_format(source), GPUMD_XYZ)
            prepared, source_format = prepare_trajectory(
                source,
                batch_size=1,
            )
            self.assertEqual(source_format, GPUMD_XYZ)
            with NetCDFReader(prepared) as trajectory:
                self.assertEqual(trajectory.read(POSITIONS).dtype, np.dtype("float32"))
                np.testing.assert_allclose(
                    trajectory.read(POSITIONS)[:, 1, 0],
                    [1.0, 1.1],
                )
                np.testing.assert_allclose(
                    trajectory.read(VELOCITIES)[0, 0],
                    [100.0, 200.0, 300.0],
                    atol=1.0e-4,
                )
                np.testing.assert_array_equal(
                    trajectory.read(TYPES)[0],
                    [0, 1],
                )

            converted = default_converted_path(source)
            first_mtime = converted.stat().st_mtime_ns
            prepare_trajectory(source)
            self.assertEqual(converted.stat().st_mtime_ns, first_mtime)

    def test_lammps_dump_sorts_atom_ids_and_converts_units(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.data"
            source.write_text(LAMMPS_TEXT, encoding="utf-8")

            self.assertEqual(detect_text_format(source), LAMMPS_DUMP)
            prepared, source_format = prepare_trajectory(
                source,
                batch_size=1,
                lammps_unit="metal",
            )
            self.assertEqual(source_format, LAMMPS_DUMP)
            with NetCDFReader(prepared) as trajectory:
                np.testing.assert_array_equal(
                    trajectory.read(ATOM_IDS),
                    [1, 2],
                )
                np.testing.assert_allclose(
                    trajectory.read(POSITIONS)[0, :, 0],
                    [1, 2],
                )
                np.testing.assert_allclose(
                    trajectory.read(VELOCITIES)[0, 0],
                    [100, 200, 300],
                )
                np.testing.assert_allclose(
                    np.diag(trajectory.read_cells()[0]),
                    [10, 20, 30],
                )

    def test_lammps_unit_change_invalidates_the_cached_conversion(self) -> None:
        with TemporaryDirectory() as directory:
            source = Path(directory) / "dump.data"
            source.write_text(LAMMPS_TEXT, encoding="utf-8")

            cached, _ = prepare_trajectory(
                source,
                compression_level=0,
                lammps_unit="metal",
            )
            with NetCDFReader(cached) as trajectory:
                metal_velocity = trajectory.read(VELOCITIES)[0, 0, 0]

            cached_again, _ = prepare_trajectory(
                source,
                compression_level=0,
                lammps_unit="real",
            )
            with NetCDFReader(cached_again) as trajectory:
                real_velocity = trajectory.read(VELOCITIES)[0, 0, 0]

        self.assertEqual(cached_again, cached)
        self.assertAlmostEqual(metal_velocity, 100.0)
        self.assertAlmostEqual(real_velocity, 100000.0)


if __name__ == "__main__":
    unittest.main()

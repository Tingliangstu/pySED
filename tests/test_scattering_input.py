from pathlib import Path

import numpy as np
import pytest

from pySED import My_Parsers
from pySED.scattering_input import ScatteringInputError, run_q_advisor_from_input


def _write_input(tmp_path, text):
    path = tmp_path / "input_SED.in"
    path.write_text(text.strip() + "\n", encoding="utf-8")
    return path


def test_input_parser_handles_utf8_bom(tmp_path):
    path = tmp_path / "input_SED.in"
    path.write_bytes(b"\xef\xbb\xbfnum_atoms = 2\n")

    params = My_Parsers.get_parse_input(str(path))

    assert params.num_atoms == 2


def test_scattering_mesh_input_writes_q_advisor_outputs(tmp_path):
    prefix = tmp_path / "mesh_report"
    input_file = _write_input(
        tmp_path,
        """
        num_atoms = 2
        total_num_steps = 16
        time_step = 1.0
        output_data_stride = 1
        prim_unitcell = 1 0 0  0 1 0  0 0 1
        supercell_dim = 4 1 1
        scattering = 1
        scattering_mode = q_advisor
        scattering_qpoints_option = mesh
        scattering_q_mesh_H = 0.20 0.30 2
        scattering_q_mesh_K = 0.0
        scattering_q_mesh_L = 0.0
        scattering_q_policy = nearest
        scattering_output_prefix = %s
        """ % prefix,
    )

    params = My_Parsers.get_parse_input(str(input_file))
    result = run_q_advisor_from_input(params)

    report = result["report"]
    assert not report.all_commensurate
    np.testing.assert_allclose(report.nearest_reduced[:, 0], [0.25, 0.25])
    assert result["efficiency"].qpoint_reduction_factor == 2.0

    for path in result["written_files"].values():
        assert Path(path).exists()

    summary = Path(result["written_files"]["summary"]).read_text(encoding="utf-8")
    assert "q-point reduction factor" in summary
    assert "suggested_diagonal_supercell = [10, 1, 1]" in summary


def test_scattering_strict_input_rejects_non_commensurate_q_but_writes_report(tmp_path):
    prefix = tmp_path / "strict_report"
    input_file = _write_input(
        tmp_path,
        """
        num_atoms = 2
        total_num_steps = 16
        time_step = 1.0
        output_data_stride = 1
        prim_unitcell = 1 0 0  0 1 0  0 0 1
        supercell_dim = 4 1 1
        scattering = 1
        scattering_mode = q_advisor
        scattering_qpoints_option = points
        scattering_qpoints = 0.20 0.0 0.0
        scattering_q_policy = strict
        scattering_output_prefix = %s
        """ % prefix,
    )

    params = My_Parsers.get_parse_input(str(input_file))
    with pytest.raises(ScatteringInputError) as exc:
        run_q_advisor_from_input(params)

    assert "non-commensurate" in str(exc.value)
    assert Path(exc.value.written_files["summary"]).exists()
    assert Path(exc.value.written_files["nearest_qpoints_reduced"]).exists()


def test_scattering_path_input_uses_q_path_steps(tmp_path):
    prefix = tmp_path / "path_report"
    input_file = _write_input(
        tmp_path,
        """
        num_atoms = 2
        total_num_steps = 16
        time_step = 1.0
        output_data_stride = 1
        prim_unitcell = 1 0 0  0 1 0  0 0 1
        supercell_dim = 4 1 1
        num_qpaths = 1
        q_path = 0.0 0.0 0.0  0.5 0.0 0.0
        scattering = 1
        scattering_mode = q_advisor
        scattering_qpoints_option = path
        scattering_q_path_steps = 2
        scattering_q_policy = strict
        scattering_output_prefix = %s
        """ % prefix,
    )

    params = My_Parsers.get_parse_input(str(input_file))
    result = run_q_advisor_from_input(params)

    report = result["report"]
    assert report.all_commensurate
    assert report.num_points == 3
    np.testing.assert_allclose(report.requested_reduced[:, 0], [0.0, 0.25, 0.5])
    np.testing.assert_allclose(
        np.loadtxt(result["written_files"]["phonopy_qpoints_reduced"]).reshape(-1, 3)[:, 0],
        [0.0, 0.25, 0.5],
    )

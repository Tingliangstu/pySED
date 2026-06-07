import csv

import numpy as np
import pytest

from pySED.eels import compute_mode_visibility_decomposition
from pySED.mode_visibility import compute_one_phonon_visibility
from pySED.visibility_diagnostics import diagnose_mode_visibility


def _two_atom_interference_case():
    primitive = np.eye(3)
    basis_positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    masses = np.ones(2)
    qpoints = np.array([[0.0, 0.0, 0.0]])
    g_vectors = np.array([[1.0, 0.0, 0.0]])
    eig = np.zeros((1, 2, 2, 3), dtype=complex)
    eig[0, 0, :, 0] = [1.0, 1.0] / np.sqrt(2.0)
    eig[0, 1, :, 0] = [1.0, -1.0] / np.sqrt(2.0)
    return primitive, basis_positions, masses, qpoints, g_vectors, eig


def test_eels_visibility_diagnostics_identify_basis_extinction():
    primitive, basis_positions, masses, qpoints, g_vectors, eig = _two_atom_interference_case()
    decomposition = compute_mode_visibility_decomposition(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
    )

    diagnostic = diagnose_mode_visibility(
        decomposition,
        atom_labels=["A", "B"],
        mode_labels=["acoustic", "optic"],
    )

    assert diagnostic.dominant_reason[0, 0, 0] == "basis_interference"
    assert diagnostic.dominant_reason[0, 0, 1] == "visible"
    assert diagnostic.basis_interference_fraction[0, 0, 0] == pytest.approx(-1.0)
    assert diagnostic.relative_visibility[0, 0, 1] == pytest.approx(1.0)

    weak = diagnostic.weak_records()
    assert len(weak) == 1
    assert weak[0]["mode_label"] == "acoustic"
    assert weak[0]["dominant_atom_label"] in ("A", "B")


def test_one_phonon_visibility_diagnostics_write_csv(tmp_path):
    primitive, basis_positions, masses, qpoints, g_vectors, eig = _two_atom_interference_case()
    visibility = compute_one_phonon_visibility(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
        experiment="custom",
    )
    diagnostic = diagnose_mode_visibility(visibility, atom_labels=["A", "B"])

    path = tmp_path / "visibility_diagnostics.csv"
    diagnostic.write_csv(path, include_all=False)

    with open(path, "r", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 1
    assert rows[0]["dominant_reason"] == "basis_interference"
    assert rows[0]["q_error_cartesian"] == "0.0"


def test_visibility_diagnostics_report_finite_size_mapping_as_reason():
    primitive, basis_positions, masses, qpoints, g_vectors, eig = _two_atom_interference_case()
    decomposition = compute_mode_visibility_decomposition(
        qpoints,
        g_vectors,
        primitive,
        basis_positions,
        masses,
        eig,
    )

    class Advice:
        error_reduced = 0.125
        error_cartesian = 0.25

    diagnostic = diagnose_mode_visibility(decomposition, q_advice=[Advice()])

    assert diagnostic.dominant_reason[0, 0, 0] == "finite_size_q_mapping"
    assert diagnostic.q_error_reduced[0] == pytest.approx(0.125)
    assert diagnostic.q_error_cartesian[0] == pytest.approx(0.25)


def test_visibility_diagnostics_reject_plain_visibility_without_decomposition():
    class PlainVisibility:
        visibility = np.zeros((1, 1, 1))

    with pytest.raises(ValueError, match="missing"):
        diagnose_mode_visibility(PlainVisibility())

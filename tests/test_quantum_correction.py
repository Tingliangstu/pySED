import numpy as np
import pytest

from pySED.quantum_correction import (
    PLANCK_MEV_THZ,
    apply_quantum_correction,
    bose_population,
    classical_to_quantum_factor,
    convert_energy_axis,
)


def test_energy_axis_conversion_thz_to_mev():
    np.testing.assert_allclose(convert_energy_axis([1.0], "thz", "meV"), [PLANCK_MEV_THZ])
    np.testing.assert_allclose(convert_energy_axis([PLANCK_MEV_THZ], "meV", "thz"), [1.0])


def test_quantum_factor_obeys_detailed_balance():
    energy = np.array([-10.0, 0.0, 10.0])
    temperature = 300.0
    factor = classical_to_quantum_factor(energy, temperature, unit="meV")

    assert factor[1] == pytest.approx(1.0)
    np.testing.assert_allclose(
        factor[0] / factor[2],
        np.exp(-10.0 / (0.08617333262145 * temperature)),
    )


def test_apply_quantum_correction_multiplies_selected_axis():
    intensity = np.ones((2, 3))
    energy = np.array([0.0, 5.0, 10.0])
    corrected = apply_quantum_correction(intensity, energy, 300.0, unit="meV", axis_index=1)
    factor = classical_to_quantum_factor(energy, 300.0, unit="meV")

    np.testing.assert_allclose(corrected[0], factor)
    np.testing.assert_allclose(corrected[1], factor)


def test_bose_population_positive_energy():
    population = bose_population(np.array([10.0]), 300.0, unit="meV")

    assert population[0] > 0

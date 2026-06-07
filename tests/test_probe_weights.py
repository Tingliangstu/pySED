import numpy as np

from pySED.probe_weights import (
    CROMER_MANN_COEFFICIENTS,
    atomic_number_weights,
    cromer_mann_form_factor,
    electron_form_factor_weights,
    mott_bethe_electron_form_factor,
    xray_form_factor_weights,
)


def test_cromer_mann_q_zero_matches_atomic_number():
    for symbol in ("C", "O", "Si", "Mo", "W"):
        q_zero = cromer_mann_form_factor(0.0, CROMER_MANN_COEFFICIENTS[symbol])
        atomic_number = atomic_number_weights([symbol])[0]
        np.testing.assert_allclose(q_zero, atomic_number, rtol=2e-3, atol=5e-3)


def test_xray_form_factor_weights_are_q_dependent():
    qpoints = np.array([[0.0, 0.0, 0.0], [8.0, 0.0, 0.0]])
    weights = xray_form_factor_weights(["C", "Si"], qpoints)

    assert weights.shape == (2, 2)
    assert np.all(weights[1] < weights[0])
    np.testing.assert_allclose(weights[0], atomic_number_weights(["C", "Si"]), atol=5e-3)


def test_xray_form_factor_missing_element_options():
    qpoints = np.array([[1.0, 0.0, 0.0]])

    zero = xray_form_factor_weights(["Unknown"], qpoints, missing="zero")
    np.testing.assert_allclose(zero, [[0.0]])

    try:
        xray_form_factor_weights(["Unknown"], qpoints, missing="raise")
    except KeyError:
        pass
    else:
        raise AssertionError("missing='raise' should reject elements without coefficients")


def test_user_supplied_cromer_mann_coefficients_override_builtin():
    qpoints = np.array([[0.0, 0.0, 0.0]])
    table = {"C": {"a": (1.0, 0.0, 0.0, 0.0), "b": (0.0, 0.0, 0.0, 0.0), "c": 0.5}}
    weights = xray_form_factor_weights(["C"], qpoints, coefficients_table=table)

    np.testing.assert_allclose(weights, [[1.5]])


def test_mott_bethe_electron_form_factor_handles_scalar_and_small_q():
    scalar = mott_bethe_electron_form_factor(0.0, 6, CROMER_MANN_COEFFICIENTS["C"])
    vector = mott_bethe_electron_form_factor(
        np.array([0.0, 2.0]),
        6,
        CROMER_MANN_COEFFICIENTS["C"],
    )

    assert isinstance(scalar, float)
    assert scalar > 0.0
    assert vector.shape == (2,)
    assert np.all(np.isfinite(vector))
    assert np.all(vector > 0.0)
    assert vector[1] != vector[0]


def test_electron_form_factor_weights_are_q_dependent():
    qpoints = np.array([[0.5, 0.0, 0.0], [4.0, 0.0, 0.0]])
    weights = electron_form_factor_weights(["C", "Si"], qpoints)

    assert weights.shape == (2, 2)
    assert np.all(np.isfinite(weights))
    assert np.all(weights > 0.0)
    assert not np.allclose(weights[0], weights[1])


def test_electron_form_factor_missing_element_options():
    qpoints = np.array([[1.0, 0.0, 0.0]])

    zero = electron_form_factor_weights(["Unknown"], qpoints, missing="zero")
    np.testing.assert_allclose(zero, [[0.0]])

    try:
        electron_form_factor_weights(["Unknown"], qpoints, missing="raise")
    except KeyError:
        pass
    else:
        raise AssertionError("missing='raise' should reject elements without electron data")

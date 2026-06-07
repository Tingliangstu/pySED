import numpy as np
import pytest

from pySED.electron_kinematics import (
    angular_resolution_to_q_sigma,
    eels_qpoints_from_angle_axis,
    electron_beta,
    electron_wavelength_angstrom,
    electron_wavenumber_per_angstrom,
    longitudinal_momentum_transfer,
    scattering_angles_to_qpoints,
)


def test_relativistic_electron_wavelength_at_200kev():
    wavelength = electron_wavelength_angstrom(200000.0)

    assert wavelength == pytest.approx(0.02508, rel=5e-4)
    assert electron_wavenumber_per_angstrom(200000.0) == pytest.approx(2.0 * np.pi / wavelength)
    assert 0.69 < electron_beta(200000.0) < 0.70


def test_scattering_angles_to_qpoints_uses_mrad_axis():
    qpoints = scattering_angles_to_qpoints(
        theta_x=np.array([0.0, 10.0]),
        theta_y=0.0,
        beam_energy_ev=200000.0,
        angle_unit="mrad",
        include_longitudinal=False,
    )

    k0 = electron_wavenumber_per_angstrom(200000.0)
    np.testing.assert_allclose(qpoints[:, 0], [0.0, k0 * 0.010])
    np.testing.assert_allclose(qpoints[:, 1:], 0.0)


def test_longitudinal_momentum_transfer_matches_hbar_v_formula():
    q_parallel = longitudinal_momentum_transfer(100.0, 200000.0)
    qpoint = scattering_angles_to_qpoints(
        theta_x=0.0,
        beam_energy_ev=200000.0,
        energy_loss_ev=100.0,
        include_longitudinal=True,
    )

    assert q_parallel == pytest.approx(0.0728, rel=5e-3)
    assert qpoint.shape == (3,)
    assert qpoint[2] == pytest.approx(q_parallel)


def test_eels_angle_axis_direction_is_normalized():
    qpoints = eels_qpoints_from_angle_axis(
        np.array([0.0, 2.0]),
        beam_energy_ev=100000.0,
        direction=(1.0, 1.0),
        angle_unit="mrad",
    )

    assert qpoints.shape == (2, 3)
    assert qpoints[1, 0] == pytest.approx(qpoints[1, 1])
    assert qpoints[1, 2] == 0.0


def test_eels_angle_axis_rejects_zero_direction():
    with pytest.raises(ValueError, match="direction"):
        eels_qpoints_from_angle_axis([0.0, 1.0], 100000.0, direction=(0.0, 0.0))


def test_angular_resolution_to_q_sigma_matches_small_angle_limit():
    k0 = electron_wavenumber_per_angstrom(200000.0)

    sigma_q = angular_resolution_to_q_sigma(2.0, 200000.0, angle_unit="mrad")

    assert sigma_q == pytest.approx(k0 * 0.002)

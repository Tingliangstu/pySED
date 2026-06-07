"""Electron-beam kinematics for calibrating q-EELS momentum axes."""

import numpy as np


ELECTRON_REST_ENERGY_EV = 510998.95000
HC_EV_ANGSTROM = 12398.419843320026
HBAR_C_EV_ANGSTROM = 1973.269804593025


def _positive_energy(value, name):
    value = float(value)
    if value <= 0.0:
        raise ValueError("%s must be positive" % name)
    return value


def electron_momentum_pc_ev(beam_energy_ev):
    """Return relativistic electron momentum times c in eV."""

    kinetic = _positive_energy(beam_energy_ev, "beam_energy_ev")
    return float(np.sqrt(kinetic * (kinetic + 2.0 * ELECTRON_REST_ENERGY_EV)))


def electron_wavelength_angstrom(beam_energy_ev):
    """Return relativistic de Broglie wavelength in Angstrom."""

    return HC_EV_ANGSTROM / electron_momentum_pc_ev(beam_energy_ev)


def electron_wavenumber_per_angstrom(beam_energy_ev):
    """Return electron wave number ``2*pi/lambda`` in inverse Angstrom."""

    return 2.0 * np.pi / electron_wavelength_angstrom(beam_energy_ev)


def electron_beta(beam_energy_ev):
    """Return relativistic electron speed as ``v/c``."""

    kinetic = _positive_energy(beam_energy_ev, "beam_energy_ev")
    pc = electron_momentum_pc_ev(kinetic)
    return pc / (kinetic + ELECTRON_REST_ENERGY_EV)


def longitudinal_momentum_transfer(energy_loss_ev, beam_energy_ev):
    """Return the small-loss longitudinal momentum transfer in inverse Angstrom.

    The approximation is ``q_parallel = DeltaE / (hbar v)`` with
    ``v = beta c``.  For vibrational losses this term is usually much smaller
    than the transverse detector momentum but it is useful for documenting the
    experiment geometry.
    """

    loss = np.asarray(energy_loss_ev, dtype=float)
    beta = electron_beta(beam_energy_ev)
    return loss / (HBAR_C_EV_ANGSTROM * beta)


def _angle_scale(unit):
    key = str(unit).lower()
    if key in ("rad", "radian", "radians"):
        return 1.0
    if key in ("mrad", "milliradian", "milliradians"):
        return 1.0e-3
    if key in ("degree", "degrees", "deg"):
        return np.pi / 180.0
    raise ValueError("angle_unit must be 'rad', 'mrad', or 'degree'")


def scattering_angles_to_qpoints(
    theta_x,
    theta_y=None,
    beam_energy_ev=200000.0,
    energy_loss_ev=0.0,
    angle_unit="rad",
    include_longitudinal=True,
):
    """Convert small electron scattering angles to Cartesian Q points.

    ``theta_x`` and ``theta_y`` are detector scattering angles.  The transverse
    momentum transfer is reported with the experimental sign convention
    ``Q_perp = k0 * theta`` where ``k0 = 2*pi/lambda``.  If
    ``include_longitudinal`` is true, the third component is
    ``DeltaE / (hbar v)``; otherwise it is zero.
    """

    theta_x = np.asarray(theta_x, dtype=float)
    if theta_y is None:
        theta_y = np.zeros_like(theta_x)
    theta_y = np.asarray(theta_y, dtype=float)
    theta_x, theta_y = np.broadcast_arrays(theta_x, theta_y)
    scale = _angle_scale(angle_unit)
    theta_x_rad = theta_x * scale
    theta_y_rad = theta_y * scale
    k0 = electron_wavenumber_per_angstrom(beam_energy_ev)

    q_z = np.zeros_like(theta_x_rad, dtype=float)
    if include_longitudinal:
        q_z = np.broadcast_to(
            longitudinal_momentum_transfer(energy_loss_ev, beam_energy_ev),
            theta_x_rad.shape,
        ).astype(float)

    return np.stack((k0 * theta_x_rad, k0 * theta_y_rad, q_z), axis=-1)


def eels_qpoints_from_angle_axis(
    angle_axis,
    beam_energy_ev,
    direction=(1.0, 0.0),
    energy_loss_ev=0.0,
    angle_unit="mrad",
    include_longitudinal=False,
):
    """Return a q-EELS line-scan Q path from a scalar scattering-angle axis."""

    axis = np.asarray(angle_axis, dtype=float)
    if axis.ndim != 1:
        raise ValueError("angle_axis must be one-dimensional")
    direction = np.asarray(direction, dtype=float)
    if direction.shape != (2,):
        raise ValueError("direction must contain two transverse components")
    norm = float(np.linalg.norm(direction))
    if norm == 0.0:
        raise ValueError("direction must be nonzero")
    direction = direction / norm
    return scattering_angles_to_qpoints(
        theta_x=axis * direction[0],
        theta_y=axis * direction[1],
        beam_energy_ev=beam_energy_ev,
        energy_loss_ev=energy_loss_ev,
        angle_unit=angle_unit,
        include_longitudinal=include_longitudinal,
    )


def angular_resolution_to_q_sigma(sigma_angle, beam_energy_ev, angle_unit="mrad"):
    """Convert an angular resolution width to transverse ``sigma_q``.

    The returned value is in inverse Angstrom and can be passed to
    :func:`pySED.compare_experiment.prepare_map_comparison` as ``sigma_q`` when
    the scattering-map q axis is a Cartesian path distance in inverse Angstrom.
    """

    sigma = abs(float(sigma_angle)) * _angle_scale(angle_unit)
    return electron_wavenumber_per_angstrom(beam_energy_ev) * sigma

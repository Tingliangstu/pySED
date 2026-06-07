"""Probe-specific scattering weights.

This module keeps scattering constants separate from the DSF estimators.  The
X-ray path uses the Cromer-Mann analytic approximation for the neutral-atom
form factor,

    f0(s) = sum_i a_i exp(-b_i s**2) + c,  s = abs(Q) / (4*pi).

``Q`` is expected in inverse Angstrom.
"""

from dataclasses import dataclass

import numpy as np


ATOMIC_NUMBERS = {
    "H": 1,
    "B": 5,
    "C": 6,
    "N": 7,
    "O": 8,
    "Na": 11,
    "Al": 13,
    "Si": 14,
    "P": 15,
    "S": 16,
    "Mo": 42,
    "Cu": 29,
    "Se": 34,
    "Zr": 40,
    "W": 74,
}


@dataclass(frozen=True)
class CromerMannCoefficients:
    """Four-Gaussian Cromer-Mann coefficients for an atomic form factor."""

    a: tuple
    b: tuple
    c: float
    source: str = "Cromer and Mann neutral-atom analytic form"

    def __post_init__(self):
        if len(self.a) != 4 or len(self.b) != 4:
            raise ValueError("Cromer-Mann coefficients require four a and four b values")


def _cm(a1, b1, a2, b2, a3, b3, a4, b4, c):
    return CromerMannCoefficients(
        a=(float(a1), float(a2), float(a3), float(a4)),
        b=(float(b1), float(b2), float(b3), float(b4)),
        c=float(c),
    )


# A compact built-in table for elements already used in pySED's neutron/X-ray
# defaults.  Users can pass their own table for ions, anomalous terms, or
# elements not included here.
CROMER_MANN_COEFFICIENTS = {
    "H": _cm(0.493002, 10.5109, 0.322912, 26.1257, 0.140191, 3.14236, 0.040810, 57.7997, 0.003038),
    "B": _cm(2.0545, 23.2185, 1.3326, 1.0210, 1.0979, 60.3498, 0.7068, 0.1403, -0.1932),
    "C": _cm(2.3100, 20.8439, 1.0200, 10.2075, 1.5886, 0.5687, 0.8650, 51.6512, 0.2156),
    "N": _cm(12.2126, 0.0057, 3.1322, 9.8933, 2.0125, 28.9975, 1.1663, 0.5826, -11.529),
    "O": _cm(3.0485, 13.2771, 2.2868, 5.7011, 1.5463, 0.3239, 0.8670, 32.9089, 0.2508),
    "Na": _cm(4.7626, 3.2850, 3.1736, 8.8422, 1.2674, 0.3136, 1.1128, 129.424, 0.6760),
    "Al": _cm(6.4202, 3.0387, 1.9002, 0.7426, 1.5936, 31.5472, 1.9646, 85.0886, 1.1151),
    "Si": _cm(6.2915, 2.4386, 3.0353, 32.3337, 1.9891, 0.6785, 1.5410, 81.6937, 1.1407),
    "P": _cm(6.4345, 1.9067, 4.1791, 27.1570, 1.7800, 0.5260, 1.4908, 68.1645, 1.1149),
    "S": _cm(6.9053, 1.4679, 5.2034, 22.2151, 1.4379, 0.2536, 1.5863, 56.1720, 0.8669),
    "Mo": _cm(3.7025, 0.2772, 17.2356, 1.0958, 12.8876, 11.0040, 3.7429, 61.6584, 4.3875),
    "Cu": _cm(13.3380, 3.5828, 7.1676, 0.2470, 5.6158, 11.3966, 1.6735, 64.8126, 1.1910),
    "Se": _cm(17.0006, 2.4098, 5.8196, 0.2726, 3.9731, 15.2372, 4.3543, 43.8163, 2.8409),
    "Zr": _cm(17.8765, 1.27618, 10.9480, 11.9160, 5.41732, 0.117622, 3.65721, 87.6627, 2.06929),
    "W": _cm(29.0818, 1.72029, 15.4300, 9.22590, 14.4327, 0.321703, 5.11982, 57.0560, 9.88750),
}


def coerce_cromer_mann_coefficients(value):
    """Return ``value`` as :class:`CromerMannCoefficients`.

    ``value`` can be a ``CromerMannCoefficients`` instance, a mapping with
    ``a``, ``b``, and ``c`` keys, or a flat nine-value iterable ordered as
    ``a1, b1, a2, b2, a3, b3, a4, b4, c``.
    """

    if isinstance(value, CromerMannCoefficients):
        return value
    if isinstance(value, dict):
        return CromerMannCoefficients(
            a=tuple(value["a"]),
            b=tuple(value["b"]),
            c=float(value["c"]),
            source=value.get("source", "user supplied"),
        )

    values = tuple(float(item) for item in value)
    if len(values) != 9:
        raise ValueError("flat Cromer-Mann coefficients must contain nine values")
    return _cm(*values)


def cromer_mann_form_factor(q_magnitude, coefficients):
    """Evaluate the Cromer-Mann form factor for ``|Q|`` in inverse Angstrom."""

    coeffs = coerce_cromer_mann_coefficients(coefficients)
    q = np.asarray(q_magnitude, dtype=float)
    s_squared = (q / (4.0 * np.pi)) ** 2
    a = np.asarray(coeffs.a, dtype=float)
    b = np.asarray(coeffs.b, dtype=float)
    values = np.sum(a.reshape((-1,) + (1,) * s_squared.ndim) * np.exp(-b.reshape((-1,) + (1,) * s_squared.ndim) * s_squared), axis=0)
    return values + coeffs.c


def mott_bethe_electron_form_factor(q_magnitude, atomic_number, coefficients, small_q=1e-8):
    """Return a relative electron scattering factor from Mott-Bethe.

    The returned value is proportional to ``(Z - f_x(Q)) / |Q|**2`` with
    ``Q`` in inverse Angstrom.  The dimensional constant is omitted because
    pySED's kinematic EELS visibility uses relative intensities by default.
    The ``Q -> 0`` limit is evaluated analytically from the Cromer-Mann
    derivative.
    """

    coeffs = coerce_cromer_mann_coefficients(coefficients)
    scalar_input = np.ndim(q_magnitude) == 0
    q = np.atleast_1d(np.asarray(q_magnitude, dtype=float))
    q_squared = q ** 2
    xray = cromer_mann_form_factor(q, coeffs)
    numerator = float(atomic_number) - xray
    factor = np.empty_like(q, dtype=float)
    small = np.abs(q) < float(small_q)

    with np.errstate(divide="ignore", invalid="ignore"):
        factor[~small] = numerator[~small] / q_squared[~small]

    a = np.asarray(coeffs.a, dtype=float)
    b = np.asarray(coeffs.b, dtype=float)
    limit = float(np.sum(a * b) / (16.0 * np.pi ** 2))
    factor[small] = limit
    if scalar_input:
        return float(factor[0])
    return factor


def atomic_number_weights(atom_types, qpoints_cartesian=None, missing="raise"):
    """Return atomic-number X-ray weights.

    This is a fallback for elements without form-factor coefficients.  It is
    useful for exploratory calculations but should not be treated as an
    experiment-quality X-ray model.
    """

    if atom_types is None:
        return None
    values = []
    for atom in atom_types:
        key = str(atom)
        if key not in ATOMIC_NUMBERS:
            if missing == "zero":
                values.append(0.0)
                continue
            raise KeyError("missing atomic number for atom type %s" % key)
        values.append(float(ATOMIC_NUMBERS[key]))

    weights = np.asarray(values, dtype=float)
    if qpoints_cartesian is None:
        return weights
    qpoints = np.atleast_2d(np.asarray(qpoints_cartesian, dtype=float))
    return np.tile(weights.reshape(1, -1), (qpoints.shape[0], 1))


def xray_form_factor_weights(
    atom_types,
    qpoints_cartesian=None,
    coefficients_table=None,
    missing="atomic_number",
):
    """Return q-dependent X-ray form-factor weights.

    Parameters
    ----------
    atom_types
        Element symbols for each atom.
    qpoints_cartesian
        ``Q`` points in inverse Angstrom.  If omitted, ``Q=0`` weights are
        returned as a one-dimensional atom array.
    coefficients_table
        Optional mapping from atom symbol to Cromer-Mann coefficients.  Entries
        override the built-in neutral-atom table.
    missing
        ``"atomic_number"`` falls back to :func:`atomic_number_weights`,
        ``"zero"`` inserts zero weights, and ``"raise"`` raises ``KeyError``.
    """

    if atom_types is None:
        return None
    if missing not in ("atomic_number", "zero", "raise"):
        raise ValueError("missing must be 'atomic_number', 'zero', or 'raise'")

    table = dict(CROMER_MANN_COEFFICIENTS)
    if coefficients_table:
        table.update(coefficients_table)

    qpoints = None
    if qpoints_cartesian is not None:
        qpoints = np.atleast_2d(np.asarray(qpoints_cartesian, dtype=float))
        if qpoints.shape[1] != 3:
            raise ValueError("qpoints_cartesian must have shape (num_qpoints, 3)")
        q_magnitude = np.linalg.norm(qpoints, axis=1)
    else:
        q_magnitude = np.asarray([0.0], dtype=float)

    output = np.zeros((q_magnitude.shape[0], len(atom_types)), dtype=float)
    for atom_index, atom in enumerate(atom_types):
        key = str(atom)
        if key in table:
            output[:, atom_index] = cromer_mann_form_factor(q_magnitude, table[key])
        elif missing == "atomic_number":
            if key not in ATOMIC_NUMBERS:
                raise KeyError("missing atomic number for atom type %s" % key)
            output[:, atom_index] = float(ATOMIC_NUMBERS[key])
        elif missing == "zero":
            output[:, atom_index] = 0.0
        else:
            raise KeyError("missing Cromer-Mann coefficients for atom type %s" % key)

    if qpoints_cartesian is None:
        return output[0]
    return output


def electron_form_factor_weights(
    atom_types,
    qpoints_cartesian,
    coefficients_table=None,
    missing="raise",
):
    """Return relative electron form-factor weights using Mott-Bethe.

    This is a neutral-atom, first-Born approximation suitable for kinematic
    EELS visibility estimates.  For quantitative electron diffraction or
    multislice calculations, users should pass a dedicated electron-scattering
    table through the higher-level EELS API.
    """

    if atom_types is None:
        return None
    if missing not in ("zero", "raise"):
        raise ValueError("missing must be 'zero' or 'raise'")

    qpoints = np.atleast_2d(np.asarray(qpoints_cartesian, dtype=float))
    if qpoints.shape[1] != 3:
        raise ValueError("qpoints_cartesian must have shape (num_qpoints, 3)")
    q_magnitude = np.linalg.norm(qpoints, axis=1)

    table = dict(CROMER_MANN_COEFFICIENTS)
    if coefficients_table:
        table.update(coefficients_table)

    output = np.zeros((q_magnitude.shape[0], len(atom_types)), dtype=float)
    for atom_index, atom in enumerate(atom_types):
        key = str(atom)
        if key in table and key in ATOMIC_NUMBERS:
            output[:, atom_index] = mott_bethe_electron_form_factor(
                q_magnitude,
                ATOMIC_NUMBERS[key],
                table[key],
            )
        elif missing == "zero":
            output[:, atom_index] = 0.0
        else:
            raise KeyError("missing electron form-factor data for atom type %s" % key)
    return output

"""Mode-level visibility diagnostics for experiment-facing scattering maps."""

from dataclasses import dataclass
import csv

import numpy as np


_DIRECTION_LABELS = ("x", "y", "z")


@dataclass
class ModeVisibilityDiagnostics:
    """Branch-level explanation of EELS/INS/IXS mode visibility.

    The arrays have shape ``(num_q, num_g, num_modes)`` unless noted
    otherwise.  The diagnostic quantities are derived from the complex
    visibility decomposition

    ``I = |sum_{b,alpha} A_{b,alpha}|^2``.
    """

    visibility: np.ndarray
    relative_visibility: np.ndarray
    atom_sum: np.ndarray
    direction_sum: np.ndarray
    basis_interference: np.ndarray
    direction_interference: np.ndarray
    basis_interference_fraction: np.ndarray
    direction_interference_fraction: np.ndarray
    dominant_atom_index: np.ndarray
    dominant_direction_index: np.ndarray
    dominant_reason: np.ndarray
    q_error_reduced: np.ndarray
    q_error_cartesian: np.ndarray
    atom_labels: list
    direction_labels: list
    mode_labels: list
    metadata: dict

    def reason_counts(self):
        """Return a dictionary counting the dominant reason labels."""

        labels, counts = np.unique(self.dominant_reason.astype(str), return_counts=True)
        return {str(label): int(count) for label, count in zip(labels, counts)}

    def to_records(self, include_all=True, relative_threshold=None):
        """Return a list of row dictionaries suitable for tables or CSV.

        Parameters
        ----------
        include_all : bool
            If ``False``, only rows with relative visibility less than
            ``relative_threshold`` are emitted.
        relative_threshold : float or None
            Relative-visibility cutoff for weak-branch rows.  When ``None``,
            the value stored in ``metadata["relative_floor"]`` is used.
        """

        cutoff = self.metadata.get("relative_floor", 1e-3)
        if relative_threshold is not None:
            cutoff = float(relative_threshold)

        records = []
        n_q, n_g, n_modes = self.visibility.shape
        for i_q in range(n_q):
            for i_g in range(n_g):
                for i_mode in range(n_modes):
                    rel = float(self.relative_visibility[i_q, i_g, i_mode])
                    if not include_all and rel >= cutoff:
                        continue
                    atom_index = int(self.dominant_atom_index[i_q, i_g, i_mode])
                    direction_index = int(self.dominant_direction_index[i_q, i_g, i_mode])
                    records.append(
                        {
                            "q_index": i_q,
                            "g_index": i_g,
                            "mode_index": i_mode,
                            "mode_label": self.mode_labels[i_mode],
                            "visibility": float(self.visibility[i_q, i_g, i_mode]),
                            "relative_visibility": rel,
                            "atom_sum": float(self.atom_sum[i_q, i_g, i_mode]),
                            "direction_sum": float(self.direction_sum[i_q, i_g, i_mode]),
                            "basis_interference": float(self.basis_interference[i_q, i_g, i_mode]),
                            "basis_interference_fraction": float(
                                self.basis_interference_fraction[i_q, i_g, i_mode]
                            ),
                            "direction_interference": float(
                                self.direction_interference[i_q, i_g, i_mode]
                            ),
                            "direction_interference_fraction": float(
                                self.direction_interference_fraction[i_q, i_g, i_mode]
                            ),
                            "dominant_atom_index": atom_index,
                            "dominant_atom_label": self.atom_labels[atom_index],
                            "dominant_direction": self.direction_labels[direction_index],
                            "dominant_reason": str(self.dominant_reason[i_q, i_g, i_mode]),
                            "q_error_reduced": float(self.q_error_reduced[i_q]),
                            "q_error_cartesian": float(self.q_error_cartesian[i_q]),
                        }
                    )
        return records

    def weak_records(self, relative_threshold=None):
        """Return only branches classified as weak by relative visibility."""

        return self.to_records(include_all=False, relative_threshold=relative_threshold)

    def write_csv(self, path, include_all=True, relative_threshold=None):
        """Write diagnostic records to a CSV table."""

        records = self.to_records(include_all=include_all, relative_threshold=relative_threshold)
        fieldnames = [
            "q_index",
            "g_index",
            "mode_index",
            "mode_label",
            "visibility",
            "relative_visibility",
            "atom_sum",
            "direction_sum",
            "basis_interference",
            "basis_interference_fraction",
            "direction_interference",
            "direction_interference_fraction",
            "dominant_atom_index",
            "dominant_atom_label",
            "dominant_direction",
            "dominant_reason",
            "q_error_reduced",
            "q_error_cartesian",
        ]
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)


def _safe_fraction(numerator, denominator):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    output = np.zeros_like(numerator, dtype=float)
    mask = np.abs(denominator) > 0.0
    output[mask] = numerator[mask] / denominator[mask]
    return output


def _relative_visibility(visibility):
    local_max = np.max(visibility, axis=2, keepdims=True)
    relative = np.zeros_like(visibility, dtype=float)
    np.divide(visibility, local_max, out=relative, where=local_max > 0.0)
    return relative


def _labels(labels, size, prefix):
    if labels is None:
        return ["%s%d" % (prefix, index) for index in range(size)]
    labels = [str(label) for label in labels]
    if len(labels) != size:
        raise ValueError("%s labels length must be %d" % (prefix, size))
    return labels


def _q_errors_from_advice(q_advice, n_q):
    q_error_reduced = np.zeros(n_q, dtype=float)
    q_error_cartesian = np.zeros(n_q, dtype=float)
    if q_advice is None:
        return q_error_reduced, q_error_cartesian
    if hasattr(q_advice, "error_reduced") and hasattr(q_advice, "error_cartesian"):
        reduced = np.asarray(q_advice.error_reduced, dtype=float)
        cartesian = np.asarray(q_advice.error_cartesian, dtype=float)
        if reduced.shape != (n_q,) or cartesian.shape != (n_q,):
            raise ValueError("q_advice error arrays must match number of q points")
        return reduced, cartesian
    if len(q_advice) != n_q:
        raise ValueError("q_advice length must match number of q points")
    for index, advice in enumerate(q_advice):
        q_error_reduced[index] = float(getattr(advice, "error_reduced", 0.0))
        q_error_cartesian[index] = float(getattr(advice, "error_cartesian", 0.0))
    return q_error_reduced, q_error_cartesian


def _classify_reasons(
    visibility,
    relative_visibility,
    atom_sum,
    direction_sum,
    basis_fraction,
    direction_fraction,
    q_error_cartesian,
    visibility_floor,
    relative_floor,
    interference_floor,
    q_error_floor,
):
    reasons = np.empty(visibility.shape, dtype=object)
    n_q, n_g, n_modes = visibility.shape
    for i_q in range(n_q):
        q_mapped = q_error_cartesian[i_q] > q_error_floor
        for i_g in range(n_g):
            for i_mode in range(n_modes):
                value = visibility[i_q, i_g, i_mode]
                relative = relative_visibility[i_q, i_g, i_mode]
                if value > visibility_floor and relative >= relative_floor:
                    reasons[i_q, i_g, i_mode] = "visible"
                    continue
                if q_mapped:
                    reasons[i_q, i_g, i_mode] = "finite_size_q_mapping"
                elif atom_sum[i_q, i_g, i_mode] <= visibility_floor:
                    reasons[i_q, i_g, i_mode] = "polarization_or_zero_weight"
                elif basis_fraction[i_q, i_g, i_mode] <= -interference_floor:
                    reasons[i_q, i_g, i_mode] = "basis_interference"
                elif direction_sum[i_q, i_g, i_mode] <= visibility_floor:
                    reasons[i_q, i_g, i_mode] = "direction_selection"
                elif direction_fraction[i_q, i_g, i_mode] <= -interference_floor:
                    reasons[i_q, i_g, i_mode] = "cartesian_interference"
                else:
                    reasons[i_q, i_g, i_mode] = "weak_visibility"
    return reasons


def diagnose_mode_visibility(
    visibility_decomposition,
    q_advice=None,
    atom_labels=None,
    mode_labels=None,
    visibility_floor=1e-14,
    relative_floor=1e-3,
    interference_floor=0.8,
    q_error_floor=1e-10,
):
    """Diagnose why each branch is visible or weak in EELS/INS/IXS.

    ``visibility_decomposition`` can be an
    :class:`pySED.eels.EELSVisibilityDecomposition` or a
    :class:`pySED.mode_visibility.OnePhononVisibility`.  It must expose
    ``visibility``, ``atom_visibility`` and ``direction_visibility``.  The
    returned object quantifies atom/basis interference and Cartesian-direction
    selection without changing the underlying scattering intensity.
    """

    required = ("visibility", "atom_visibility", "direction_visibility")
    missing = [name for name in required if not hasattr(visibility_decomposition, name)]
    if missing:
        raise ValueError("visibility decomposition is missing: %s" % ", ".join(missing))

    visibility = np.asarray(visibility_decomposition.visibility, dtype=float)
    atom_visibility = np.asarray(visibility_decomposition.atom_visibility, dtype=float)
    direction_visibility = np.asarray(visibility_decomposition.direction_visibility, dtype=float)
    if visibility.ndim != 3:
        raise ValueError("visibility must have shape (num_q, num_g, num_modes)")
    if atom_visibility.shape[:3] != visibility.shape:
        raise ValueError("atom_visibility leading dimensions must match visibility")
    if direction_visibility.shape != visibility.shape + (3,):
        raise ValueError("direction_visibility must have shape (num_q, num_g, num_modes, 3)")

    atom_sum = np.sum(atom_visibility, axis=3)
    direction_sum = np.sum(direction_visibility, axis=3)
    basis_interference = visibility - atom_sum
    direction_interference = visibility - direction_sum
    basis_fraction = _safe_fraction(basis_interference, atom_sum)
    direction_fraction = _safe_fraction(direction_interference, direction_sum)
    relative = _relative_visibility(visibility)
    dominant_atom = np.argmax(atom_visibility, axis=3)
    dominant_direction = np.argmax(direction_visibility, axis=3)
    q_error_reduced, q_error_cartesian = _q_errors_from_advice(q_advice, visibility.shape[0])

    reasons = _classify_reasons(
        visibility,
        relative,
        atom_sum,
        direction_sum,
        basis_fraction,
        direction_fraction,
        q_error_cartesian,
        float(visibility_floor),
        float(relative_floor),
        float(interference_floor),
        float(q_error_floor),
    )

    n_modes = visibility.shape[2]
    n_atoms = atom_visibility.shape[3]
    metadata = {
        "source_model": getattr(visibility_decomposition, "metadata", {}).get("model"),
        "visibility_floor": float(visibility_floor),
        "relative_floor": float(relative_floor),
        "interference_floor": float(interference_floor),
        "q_error_floor": float(q_error_floor),
        "basis_interference": "I_total - sum_b |sum_alpha A_balpha|^2",
        "direction_interference": "I_total - sum_alpha |sum_b A_balpha|^2",
    }

    return ModeVisibilityDiagnostics(
        visibility=visibility,
        relative_visibility=relative,
        atom_sum=atom_sum,
        direction_sum=direction_sum,
        basis_interference=basis_interference,
        direction_interference=direction_interference,
        basis_interference_fraction=basis_fraction,
        direction_interference_fraction=direction_fraction,
        dominant_atom_index=dominant_atom,
        dominant_direction_index=dominant_direction,
        dominant_reason=reasons,
        q_error_reduced=q_error_reduced,
        q_error_cartesian=q_error_cartesian,
        atom_labels=_labels(atom_labels, n_atoms, "atom"),
        direction_labels=list(_DIRECTION_LABELS),
        mode_labels=_labels(mode_labels, n_modes, "mode"),
        metadata=metadata,
    )

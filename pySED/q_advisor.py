"""Q-point utilities for finite periodic MD trajectories.

The routines in this module make the supercell commensurability condition
explicit.  For a primitive cell ``p`` and a simulation supercell ``S`` with
``S = P @ p``, a reduced q point is compatible with the trajectory if

    q_red @ P.T is an integer vector.

This is the same condition used by the existing Brillouin-zone path builder,
but exposed as reusable API for SED, DSF, EELS, and experiment comparison.
"""

from dataclasses import dataclass
from fractions import Fraction
import csv
import itertools
import math

import numpy as np


def _lcm(a, b):
    return abs(a * b) // math.gcd(a, b)


@dataclass(frozen=True)
class QPointAdvice:
    """Result of mapping one requested q point to the commensurate grid."""

    requested_reduced: np.ndarray
    nearest_reduced: np.ndarray
    requested_cartesian: np.ndarray
    nearest_cartesian: np.ndarray
    is_commensurate: bool
    error_reduced: float
    error_cartesian: float
    integer_supercell_index: np.ndarray


@dataclass(frozen=True)
class QPointMesh:
    """Selected-Q mesh in reduced reciprocal coordinates."""

    qpoints_reduced: np.ndarray
    axes: tuple
    mesh_shape: tuple

    def reshape_values(self, values):
        """Reshape per-Q values back onto the selected-Q mesh."""

        arr = np.asarray(values)
        if arr.shape[0] != self.qpoints_reduced.shape[0]:
            raise ValueError("values first dimension must match number of q points")
        return arr.reshape(tuple(self.mesh_shape) + arr.shape[1:])


@dataclass(frozen=True)
class QPathSegments:
    """Selected-Q path assembled from one or more linear segments."""

    qpoints_reduced: np.ndarray
    path_axis: np.ndarray
    segment_start_indices: np.ndarray
    segment_end_indices: np.ndarray


@dataclass(frozen=True)
class SelectedQEfficiencyReport:
    """Cost estimate for selected-Q evaluation relative to a full q grid."""

    num_selected_qpoints: int
    num_full_qpoints: int
    selected_fraction: float
    qpoint_reduction_factor: float
    num_frames: int = None
    num_atoms: int = None
    selected_phase_evaluations: float = None
    full_phase_evaluations: float = None
    selected_fft_work_units: float = None
    full_fft_work_units: float = None

    def summary(self):
        """Return a JSON-serializable cost summary."""

        return {
            "num_selected_qpoints": self.num_selected_qpoints,
            "num_full_qpoints": self.num_full_qpoints,
            "selected_fraction": self.selected_fraction,
            "qpoint_reduction_factor": self.qpoint_reduction_factor,
            "num_frames": self.num_frames,
            "num_atoms": self.num_atoms,
            "selected_phase_evaluations": self.selected_phase_evaluations,
            "full_phase_evaluations": self.full_phase_evaluations,
            "selected_fft_work_units": self.selected_fft_work_units,
            "full_fft_work_units": self.full_fft_work_units,
        }

    def to_table(self, precision=6):
        """Return a compact plain-text table for logs or methods notes."""

        fmt = "%%.%df" % int(precision)

        def val(item):
            if item is None:
                return "n/a"
            if isinstance(item, (int, np.integer)):
                return str(int(item))
            return fmt % float(item)

        rows = [
            ("selected q points", val(self.num_selected_qpoints)),
            ("full commensurate grid points", val(self.num_full_qpoints)),
            ("selected fraction", val(self.selected_fraction)),
            ("q-point reduction factor", val(self.qpoint_reduction_factor)),
            ("selected phase evaluations", val(self.selected_phase_evaluations)),
            ("full-grid phase evaluations", val(self.full_phase_evaluations)),
            ("selected FFT work units", val(self.selected_fft_work_units)),
            ("full-grid FFT work units", val(self.full_fft_work_units)),
        ]
        width = max(len(label) for label, _ in rows)
        return "\n".join("%s  %s" % (label.ljust(width), value) for label, value in rows)


@dataclass(frozen=True)
class QAdvisorReport:
    """Vectorized report for an experimental q or Q path.

    The report keeps the requested points, the nearest commensurate points, the
    finite-supercell integer indices, and a diagonal-supercell recommendation in
    one object so scripts can print, export, or archive the q-selection decision
    before running SED, DSF, or EELS calculations.
    """

    coordinate_system: str
    preserve_image: bool
    requested_reduced: np.ndarray
    nearest_reduced: np.ndarray
    requested_cartesian: np.ndarray
    nearest_cartesian: np.ndarray
    is_commensurate: np.ndarray
    error_reduced: np.ndarray
    error_cartesian: np.ndarray
    integer_supercell_indices: np.ndarray
    suggested_diagonal_supercell: np.ndarray

    @property
    def num_points(self):
        return int(self.requested_reduced.shape[0])

    @property
    def num_commensurate(self):
        return int(np.count_nonzero(self.is_commensurate))

    @property
    def all_commensurate(self):
        return bool(np.all(self.is_commensurate))

    @property
    def max_error_reduced(self):
        return float(np.max(self.error_reduced)) if self.error_reduced.size else 0.0

    @property
    def max_error_cartesian(self):
        return float(np.max(self.error_cartesian)) if self.error_cartesian.size else 0.0

    def summary(self):
        """Return scalar report metadata useful for logs or JSON exports."""

        return {
            "coordinate_system": self.coordinate_system,
            "preserve_image": self.preserve_image,
            "num_points": self.num_points,
            "num_commensurate": self.num_commensurate,
            "all_commensurate": self.all_commensurate,
            "max_error_reduced": self.max_error_reduced,
            "max_error_cartesian": self.max_error_cartesian,
            "suggested_diagonal_supercell": self.suggested_diagonal_supercell.tolist(),
        }

    def to_records(self, flat=False):
        """Return one serializable record per requested point."""

        records = []
        for index in range(self.num_points):
            record = {
                "index": int(index),
                "is_commensurate": bool(self.is_commensurate[index]),
                "requested_reduced": self.requested_reduced[index].tolist(),
                "nearest_reduced": self.nearest_reduced[index].tolist(),
                "requested_cartesian": self.requested_cartesian[index].tolist(),
                "nearest_cartesian": self.nearest_cartesian[index].tolist(),
                "error_reduced": float(self.error_reduced[index]),
                "error_cartesian": float(self.error_cartesian[index]),
                "integer_supercell_index": self.integer_supercell_indices[index].astype(int).tolist(),
            }
            if flat:
                record = _flatten_report_record(record)
            records.append(record)
        return records

    def to_dict(self):
        """Return a JSON-serializable report dictionary."""

        return {"summary": self.summary(), "points": self.to_records()}

    def write_csv(self, path):
        """Write the q-advisor point table to ``path`` as CSV."""

        records = self.to_records(flat=True)
        fieldnames = list(records[0].keys()) if records else [
            "index",
            "is_commensurate",
            "error_reduced",
            "error_cartesian",
        ]
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)

    def to_table(self, precision=6):
        """Return a compact plain-text table for terminal logs."""

        fmt = "%%.%df" % int(precision)

        def vec(values):
            return "[" + ", ".join(fmt % float(x) for x in values) + "]"

        rows = [
            [
                str(index),
                str(bool(self.is_commensurate[index])),
                vec(self.requested_reduced[index]),
                vec(self.nearest_reduced[index]),
                fmt % float(self.error_reduced[index]),
                fmt % float(self.error_cartesian[index]),
                "[" + ", ".join(str(int(x)) for x in self.integer_supercell_indices[index]) + "]",
            ]
            for index in range(self.num_points)
        ]
        headers = ["idx", "ok", "requested_red", "nearest_red", "err_red", "err_cart", "supercell_n"]
        widths = [
            max(len(headers[col]), *(len(row[col]) for row in rows)) if rows else len(headers[col])
            for col in range(len(headers))
        ]
        lines = [
            "  ".join(headers[col].ljust(widths[col]) for col in range(len(headers))),
            "  ".join("-" * width for width in widths),
        ]
        lines.extend(
            "  ".join(row[col].ljust(widths[col]) for col in range(len(headers)))
            for row in rows
        )
        lines.append("suggested_diagonal_supercell = %s" % self.suggested_diagonal_supercell.tolist())
        return "\n".join(lines)


@dataclass(frozen=True)
class ExtendedZoneQPathPlan:
    """Pairwise ``Q = q + G`` plan for extended-zone scattering paths."""

    report: QAdvisorReport
    Q_reduced: np.ndarray
    Q_cartesian: np.ndarray
    qpoints_reduced: np.ndarray
    qpoints_cartesian: np.ndarray
    g_vectors_reduced: np.ndarray
    g_vectors_cartesian: np.ndarray

    @property
    def num_points(self):
        return int(self.Q_reduced.shape[0])

    def to_records(self, flat=False):
        records = []
        report_records = self.report.to_records(flat=False)
        for index, record in enumerate(report_records):
            enriched = dict(record)
            enriched.update(
                {
                    "Q_reduced": self.Q_reduced[index].tolist(),
                    "Q_cartesian": self.Q_cartesian[index].tolist(),
                    "folded_q_reduced": self.qpoints_reduced[index].tolist(),
                    "folded_q_cartesian": self.qpoints_cartesian[index].tolist(),
                    "G_reduced": self.g_vectors_reduced[index].astype(int).tolist(),
                    "G_cartesian": self.g_vectors_cartesian[index].tolist(),
                }
            )
            if flat:
                enriched = _flatten_extended_zone_record(enriched)
            records.append(enriched)
        return records

    def to_dict(self):
        return {"summary": self.report.summary(), "points": self.to_records()}

    def write_phonopy_qpoints(self, path):
        """Write folded commensurate q points for phonopy eigenvectors."""

        np.savetxt(path, self.qpoints_reduced, fmt="%.10f")

    def write_csv(self, path):
        """Write the extended-zone path plan to ``path`` as CSV."""

        records = self.to_records(flat=True)
        fieldnames = list(records[0].keys()) if records else ["index"]
        with open(path, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records)


def _flatten_report_record(record):
    flat = {
        "index": record["index"],
        "is_commensurate": record["is_commensurate"],
        "error_reduced": record["error_reduced"],
        "error_cartesian": record["error_cartesian"],
    }
    for key in ("requested_reduced", "nearest_reduced", "requested_cartesian", "nearest_cartesian"):
        for axis, value in enumerate(record[key]):
            flat["%s_%d" % (key, axis + 1)] = value
    for axis, value in enumerate(record["integer_supercell_index"]):
        flat["integer_supercell_index_%d" % (axis + 1)] = value
    return flat


def _flatten_extended_zone_record(record):
    flat = _flatten_report_record(record)
    for key in ("Q_reduced", "Q_cartesian", "folded_q_reduced", "folded_q_cartesian", "G_cartesian"):
        for axis, value in enumerate(record[key]):
            flat["%s_%d" % (key, axis + 1)] = value
    for axis, value in enumerate(record["G_reduced"]):
        flat["G_reduced_%d" % (axis + 1)] = value
    return flat


def compute_repetition_matrix(primitive_cell, supercell_cell, atol=1e-4):
    """Return integer ``P`` such that ``supercell_cell = P @ primitive_cell``."""

    primitive = np.asarray(primitive_cell, dtype=float)
    supercell = np.asarray(supercell_cell, dtype=float)
    if primitive.shape != (3, 3) or supercell.shape != (3, 3):
        raise ValueError("primitive_cell and supercell_cell must be 3x3 matrices")

    p_float = np.linalg.solve(primitive.T, supercell.T).T
    p_int = np.rint(p_float).astype(int)
    if not np.allclose(p_float, p_int, atol=atol):
        raise ValueError("supercell is not an integer multiple of primitive_cell")
    return p_int


def reciprocal_lattice(cell):
    """Return reciprocal lattice vectors as rows, including the 2*pi factor."""

    cell = np.asarray(cell, dtype=float)
    if cell.shape != (3, 3):
        raise ValueError("cell must be a 3x3 matrix")
    return 2.0 * np.pi * np.linalg.inv(cell.T)


def wrap_reduced(q_red, centered=False):
    """Wrap reduced coordinates into [0, 1) or [-0.5, 0.5)."""

    q = np.asarray(q_red, dtype=float)
    wrapped = q - np.floor(q)
    if centered:
        wrapped = wrapped - np.round(wrapped)
    return np.where(np.isclose(wrapped, 0.0, atol=1e-14), 0.0, wrapped)


def reduced_delta(q1, q2):
    """Shortest periodic difference between reduced reciprocal coordinates."""

    dq = np.asarray(q1, dtype=float) - np.asarray(q2, dtype=float)
    return dq - np.round(dq)


def split_extended_zone_qpoints(qpoints_reduced):
    """Split reduced extended-zone ``Q`` points into folded ``q`` and ``G``.

    The folded q points use the centered first Brillouin-zone image.  The
    returned G vectors are integer reduced reciprocal-lattice vectors satisfying
    ``Q = q + G`` within floating-point tolerance.
    """

    Q = np.atleast_2d(np.asarray(qpoints_reduced, dtype=float))
    if Q.shape[1] != 3:
        raise ValueError("qpoints_reduced must have shape (num_qpoints, 3)")
    q_folded = wrap_reduced(Q, centered=True)
    g_vectors = Q - q_folded
    g_integer = np.rint(g_vectors).astype(int)
    if not np.allclose(g_vectors, g_integer, atol=1e-8):
        raise ValueError("extended-zone split did not produce integer G vectors")
    return q_folded, g_integer


def qpoints_from_path_axis(q_axis, start_qpoint, end_qpoint):
    """Interpolate vector Q points along a calibrated scalar q-path axis.

    Experimental maps often store only a scalar horizontal axis, for example a
    path distance or detector coordinate after calibration.  This helper maps
    each scalar axis value linearly between two vector endpoints.  The returned
    vectors use the same coordinate system as ``start_qpoint`` and
    ``end_qpoint``.
    """

    axis = np.asarray(q_axis, dtype=float)
    start = np.asarray(start_qpoint, dtype=float)
    end = np.asarray(end_qpoint, dtype=float)
    if axis.ndim != 1:
        raise ValueError("q_axis must be one-dimensional")
    if axis.size < 1:
        raise ValueError("q_axis must contain at least one point")
    if start.shape != (3,) or end.shape != (3,):
        raise ValueError("start_qpoint and end_qpoint must have shape (3,)")
    if axis.size == 1:
        return start.reshape(1, 3)

    diffs = np.diff(axis)
    if np.any(diffs == 0.0) or not (np.all(diffs > 0.0) or np.all(diffs < 0.0)):
        raise ValueError("q_axis must be strictly monotonic")
    span = axis[-1] - axis[0]
    fractions = (axis - axis[0]) / span
    return start.reshape(1, 3) + fractions.reshape(-1, 1) * (end - start).reshape(1, 3)


def _mesh_axis(spec):
    arr = np.asarray(spec, dtype=float)
    if arr.ndim == 0:
        return np.asarray([float(arr)], dtype=float)
    arr = arr.ravel()
    if arr.size == 1:
        return np.asarray([float(arr[0])], dtype=float)
    if arr.size != 3:
        raise ValueError("mesh axis spec must be a scalar or [start, stop, num]")
    num = int(round(arr[2]))
    if num < 1 or not np.isclose(arr[2], num):
        raise ValueError("mesh axis num must be a positive integer")
    if num == 1:
        return np.asarray([float(arr[0])], dtype=float)
    return np.linspace(float(arr[0]), float(arr[1]), num)


def qpoints_from_mesh(h_axis, k_axis, l_axis):
    """Generate a pynamic-style selected-Q mesh in reduced coordinates.

    Each axis may be a scalar or ``[start, stop, num]``.  The function does not
    enforce commensurability; pass the returned ``qpoints_reduced`` to
    :class:`QAdvisor` before using the points with finite periodic
    trajectories.
    """

    h = _mesh_axis(h_axis)
    k = _mesh_axis(k_axis)
    l = _mesh_axis(l_axis)
    H, K, L = np.meshgrid(h, k, l, indexing="ij")
    qpoints = np.column_stack([H.ravel(), K.ravel(), L.ravel()])
    return QPointMesh(
        qpoints_reduced=qpoints,
        axes=(h, k, l),
        mesh_shape=(h.size, k.size, l.size),
    )


def qpoints_from_path_segments(starts, ends, steps, endpoint=False):
    """Generate selected Q points along one or more reduced-coordinate paths.

    This follows the selected-path convention used by pynamic-structure-factor:
    each segment is sampled with ``steps`` points; by default the segment end
    point is omitted to avoid duplicate vertices when consecutive segments
    share endpoints.  Use ``endpoint=True`` when every segment should include
    its final point.
    """

    starts = np.atleast_2d(np.asarray(starts, dtype=float))
    ends = np.atleast_2d(np.asarray(ends, dtype=float))
    steps = np.asarray(steps, dtype=int).ravel()
    if starts.shape != ends.shape or starts.shape[1] != 3:
        raise ValueError("starts and ends must both have shape (num_segments, 3)")
    if steps.shape != (starts.shape[0],):
        raise ValueError("steps length must match number of path segments")
    if np.any(steps < 1):
        raise ValueError("all path segment steps must be positive")

    qpoints = []
    path_axis = []
    start_indices = []
    end_indices = []
    distance = 0.0
    previous_end = None
    for seg_index, (start, end, count) in enumerate(zip(starts, ends, steps)):
        start_indices.append(len(qpoints))
        sample_count = int(count)
        fractions = np.linspace(0.0, 1.0, sample_count if endpoint else sample_count + 1)
        if not endpoint:
            fractions = fractions[:-1]
        segment = start.reshape(1, 3) + fractions.reshape(-1, 1) * (end - start).reshape(1, 3)
        if previous_end is not None:
            distance += float(np.linalg.norm(start - previous_end))
        segment_length = float(np.linalg.norm(end - start))
        segment_axis = distance + fractions * segment_length
        qpoints.extend(segment)
        path_axis.extend(segment_axis)
        if len(qpoints) == start_indices[-1]:
            end_indices.append(start_indices[-1])
        else:
            end_indices.append(len(qpoints) - 1)
        distance += segment_length
        previous_end = end

    return QPathSegments(
        qpoints_reduced=np.asarray(qpoints, dtype=float),
        path_axis=np.asarray(path_axis, dtype=float),
        segment_start_indices=np.asarray(start_indices, dtype=int),
        segment_end_indices=np.asarray(end_indices, dtype=int),
    )


def estimate_selected_q_efficiency(
    num_selected_qpoints,
    num_full_qpoints,
    num_frames=None,
    num_atoms=None,
):
    """Estimate selected-Q savings relative to evaluating a full q grid.

    The estimate follows the direct DSF cost model
    ``O(N_Q * N_t * N) + O(N_Q * N_t * log(N_t))``.  The phase-evaluation
    counts are exact under this model when ``num_frames`` and ``num_atoms`` are
    supplied; FFT work units are proportional estimates.
    """

    selected = int(num_selected_qpoints)
    full = int(num_full_qpoints)
    if selected < 1:
        raise ValueError("num_selected_qpoints must be positive")
    if full < 1:
        raise ValueError("num_full_qpoints must be positive")

    selected_fraction = float(selected / full)
    reduction = float(full / selected)
    frames = None if num_frames is None else int(num_frames)
    atoms = None if num_atoms is None else int(num_atoms)
    if frames is not None and frames < 1:
        raise ValueError("num_frames must be positive")
    if atoms is not None and atoms < 1:
        raise ValueError("num_atoms must be positive")

    selected_phase = None
    full_phase = None
    if frames is not None and atoms is not None:
        selected_phase = float(selected * frames * atoms)
        full_phase = float(full * frames * atoms)

    selected_fft = None
    full_fft = None
    if frames is not None:
        fft_units_per_q = float(frames * math.log2(max(frames, 2)))
        selected_fft = float(selected * fft_units_per_q)
        full_fft = float(full * fft_units_per_q)

    return SelectedQEfficiencyReport(
        num_selected_qpoints=selected,
        num_full_qpoints=full,
        selected_fraction=selected_fraction,
        qpoint_reduction_factor=reduction,
        num_frames=frames,
        num_atoms=atoms,
        selected_phase_evaluations=selected_phase,
        full_phase_evaluations=full_phase,
        selected_fft_work_units=selected_fft,
        full_fft_work_units=full_fft,
    )


class QAdvisor:
    """Construct, validate, and map q points for a finite MD supercell."""

    def __init__(self, primitive_cell, supercell_cell, atol=1e-8):
        self.primitive_cell = np.asarray(primitive_cell, dtype=float)
        self.supercell_cell = np.asarray(supercell_cell, dtype=float)
        self.atol = atol
        self.P = compute_repetition_matrix(
            self.primitive_cell, self.supercell_cell, atol=max(atol, 1e-4)
        )
        self.reciprocal_primitive = reciprocal_lattice(self.primitive_cell)
        self.num_commensurate_qpoints = int(abs(round(np.linalg.det(self.P))))

    def reduced_to_cartesian(self, qpoints_reduced):
        q = np.asarray(qpoints_reduced, dtype=float)
        return np.atleast_2d(q) @ self.reciprocal_primitive

    def cartesian_to_reduced(self, qpoints_cartesian):
        q = np.asarray(qpoints_cartesian, dtype=float)
        return np.linalg.solve(self.reciprocal_primitive.T, np.atleast_2d(q).T).T

    def supercell_indices(self, qpoints_reduced):
        q = np.atleast_2d(np.asarray(qpoints_reduced, dtype=float))
        return q @ self.P.T

    def is_commensurate(self, qpoint_reduced, atol=None):
        indices = self.supercell_indices(qpoint_reduced)[0]
        return bool(np.allclose(indices, np.rint(indices), atol=atol or self.atol))

    def validate_qpoints(self, qpoints_reduced, atol=None):
        q = np.atleast_2d(np.asarray(qpoints_reduced, dtype=float))
        indices = self.supercell_indices(q)
        return np.allclose(indices, np.rint(indices), atol=atol or self.atol, rtol=0.0)

    def estimate_selected_q_efficiency(
        self,
        qpoints_reduced=None,
        num_selected_qpoints=None,
        num_frames=None,
        num_atoms=None,
    ):
        """Estimate selected-Q savings against this supercell's full q grid."""

        if qpoints_reduced is None and num_selected_qpoints is None:
            raise ValueError("supply qpoints_reduced or num_selected_qpoints")
        if qpoints_reduced is not None:
            selected = np.atleast_2d(np.asarray(qpoints_reduced, dtype=float)).shape[0]
        else:
            selected = int(num_selected_qpoints)
        return estimate_selected_q_efficiency(
            selected,
            self.num_commensurate_qpoints,
            num_frames=num_frames,
            num_atoms=num_atoms,
        )

    def commensurate_grid(self, centered=False, max_bound=256):
        """Return all distinct commensurate q points modulo primitive G vectors.

        For diagonal supercells this is the usual ``i/Nx, j/Ny, k/Nz`` grid.  For
        non-diagonal integer repetition matrices, a bounded integer enumeration is
        used and stopped once ``abs(det(P))`` unique q points are found.
        """

        target = self.num_commensurate_qpoints
        if target < 1:
            raise ValueError("invalid repetition matrix with zero determinant")

        inv_pt = np.linalg.inv(self.P.T)
        base_bound = max(1, int(np.max(np.abs(self.P))))
        seen = {}
        for bound in range(base_bound, max_bound + 1):
            ranges = [range(-bound, bound + 1)] * 3
            for m in itertools.product(*ranges):
                q = np.asarray(m, dtype=float) @ inv_pt
                q = wrap_reduced(q, centered=centered)
                key = tuple(np.round(wrap_reduced(q, centered=False), 12))
                if key not in seen:
                    seen[key] = q
                    if len(seen) == target:
                        arr = np.array(list(seen.values()))
                        return arr[np.lexsort((arr[:, 2], arr[:, 1], arr[:, 0]))]
        raise RuntimeError("failed to enumerate all commensurate q points")

    def nearest_commensurate(self, qpoint_reduced, grid=None, preserve_image=False):
        q = np.asarray(qpoint_reduced, dtype=float)
        if grid is None:
            grid = self.commensurate_grid(centered=False)

        q_search = wrap_reduced(q, centered=False) if preserve_image else q
        deltas = np.array([reduced_delta(q_search, candidate) for candidate in grid])
        cart_deltas = deltas @ self.reciprocal_primitive
        distances = np.linalg.norm(cart_deltas, axis=1)
        idx = int(np.argmin(distances))
        nearest = q - deltas[idx] if preserve_image else grid[idx]

        requested_cart = self.reduced_to_cartesian(q)[0]
        nearest_cart = self.reduced_to_cartesian(nearest)[0]
        supercell_index = self.supercell_indices(nearest)[0]
        return QPointAdvice(
            requested_reduced=np.array(q, dtype=float),
            nearest_reduced=np.array(nearest, dtype=float),
            requested_cartesian=requested_cart,
            nearest_cartesian=nearest_cart,
            is_commensurate=self.is_commensurate(q),
            error_reduced=float(np.linalg.norm(reduced_delta(q, nearest))),
            error_cartesian=float(distances[idx]),
            integer_supercell_index=np.rint(supercell_index).astype(int),
        )

    def advise_qpoints(self, qpoints_reduced, preserve_image=False):
        grid = self.commensurate_grid(centered=False)
        return [
            self.nearest_commensurate(q, grid=grid, preserve_image=preserve_image)
            for q in np.atleast_2d(qpoints_reduced)
        ]

    def nearest_commensurate_cartesian(self, qpoint_cartesian, grid=None, preserve_image=False):
        """Map one Cartesian experimental Q point to the nearest allowed q point."""

        q_reduced = self.cartesian_to_reduced(qpoint_cartesian)[0]
        return self.nearest_commensurate(q_reduced, grid=grid, preserve_image=preserve_image)

    def advise_cartesian_qpoints(self, qpoints_cartesian, preserve_image=False):
        """Map Cartesian experimental Q points to nearest commensurate q points."""

        q_reduced = self.cartesian_to_reduced(qpoints_cartesian)
        return self.advise_qpoints(q_reduced, preserve_image=preserve_image)

    def advise_experimental_path(
        self,
        qpoints,
        coordinates="reduced",
        preserve_image=True,
        suggest_supercell=True,
    ):
        """Return a full commensurability report for an experimental q path.

        Parameters
        ----------
        qpoints
            Requested q or Q points, either reduced coordinates in the primitive
            reciprocal basis or Cartesian vectors in inverse Angstrom.
        coordinates
            ``"reduced"``/``"fractional"``/``"crystal"`` or
            ``"cartesian"``/``"angstrom^-1"``.
        preserve_image
            Preserve the extended-zone reciprocal-lattice image when mapping to
            the nearest commensurate point.  This should normally stay ``True``
            for EELS, INS, and IXS comparison because form factors and basis
            interference depend on the full ``Q = q + G``.
        suggest_supercell
            If ``True``, recommend diagonal primitive-cell repeats that would
            make all requested reduced coordinates commensurate.
        """

        coordinate_key = str(coordinates).lower()
        if coordinate_key in ("reduced", "fractional", "crystal"):
            requested_reduced = np.atleast_2d(np.asarray(qpoints, dtype=float))
            if requested_reduced.shape[1] != 3:
                raise ValueError("reduced qpoints must have shape (num_qpoints, 3)")
            coordinate_system = "reduced"
        elif coordinate_key in ("cartesian", "cart", "angstrom^-1", "1/angstrom"):
            requested_cartesian = np.atleast_2d(np.asarray(qpoints, dtype=float))
            if requested_cartesian.shape[1] != 3:
                raise ValueError("Cartesian qpoints must have shape (num_qpoints, 3)")
            requested_reduced = self.cartesian_to_reduced(requested_cartesian)
            coordinate_system = "cartesian"
        else:
            raise ValueError("coordinates must be 'reduced' or 'cartesian'")

        advice = self.advise_qpoints(requested_reduced, preserve_image=preserve_image)
        suggested = (
            suggest_diagonal_supercell(requested_reduced)
            if suggest_supercell
            else np.ones(3, dtype=int)
        )
        return QAdvisorReport(
            coordinate_system=coordinate_system,
            preserve_image=bool(preserve_image),
            requested_reduced=np.asarray([item.requested_reduced for item in advice], dtype=float),
            nearest_reduced=np.asarray([item.nearest_reduced for item in advice], dtype=float),
            requested_cartesian=np.asarray([item.requested_cartesian for item in advice], dtype=float),
            nearest_cartesian=np.asarray([item.nearest_cartesian for item in advice], dtype=float),
            is_commensurate=np.asarray([item.is_commensurate for item in advice], dtype=bool),
            error_reduced=np.asarray([item.error_reduced for item in advice], dtype=float),
            error_cartesian=np.asarray([item.error_cartesian for item in advice], dtype=float),
            integer_supercell_indices=np.asarray(
                [item.integer_supercell_index for item in advice],
                dtype=int,
            ),
            suggested_diagonal_supercell=np.asarray(suggested, dtype=int),
        )

    def plan_extended_zone_q_path(
        self,
        qpoints,
        coordinates="reduced",
        q_policy="strict",
        max_error_reduced=None,
        max_error_cartesian=None,
    ):
        """Plan a pairwise extended-zone path for EELS, INS, or IXS.

        The returned :class:`ExtendedZoneQPathPlan` contains the nearest allowed
        extended-zone :math:`Q_i`, the folded commensurate :math:`q_i`, and the
        integer reciprocal vector :math:`G_i` for each experimental point:

        ``Q_i = q_i + G_i``.
        """

        report = self.advise_experimental_path(
            qpoints,
            coordinates=coordinates,
            preserve_image=True,
        )
        policy = str(q_policy).lower()
        if policy not in ("strict", "nearest"):
            raise ValueError("q_policy must be 'strict' or 'nearest'")
        if policy == "strict" and not report.all_commensurate:
            raise ValueError("experimental Q path contains non-commensurate points")
        if policy == "nearest":
            if max_error_reduced is not None and report.max_error_reduced > max_error_reduced:
                raise ValueError("nearest commensurate q exceeds max_error_reduced")
            if max_error_cartesian is not None and report.max_error_cartesian > max_error_cartesian:
                raise ValueError("nearest commensurate q exceeds max_error_cartesian")

        Q_reduced = report.nearest_reduced
        q_folded, g_vectors = split_extended_zone_qpoints(Q_reduced)
        return ExtendedZoneQPathPlan(
            report=report,
            Q_reduced=Q_reduced,
            Q_cartesian=report.nearest_cartesian,
            qpoints_reduced=q_folded,
            qpoints_cartesian=self.reduced_to_cartesian(q_folded),
            g_vectors_reduced=g_vectors,
            g_vectors_cartesian=self.reduced_to_cartesian(g_vectors),
        )

    def advise_map_q_axis(
        self,
        scattering_map,
        start_qpoint,
        end_qpoint,
        coordinates="cartesian",
        preserve_image=True,
        suggest_supercell=True,
    ):
        """Advise Q points implied by a calibrated experimental map q axis.

        ``scattering_map`` may be any object with a one-dimensional ``q_axis``
        attribute, including :class:`pySED.compare_experiment.ScatteringMap`.
        The vector endpoints define the reciprocal-space path corresponding to
        the first and last q-axis values.
        """

        if not hasattr(scattering_map, "q_axis"):
            raise TypeError("scattering_map must expose a q_axis attribute")
        qpoints = qpoints_from_path_axis(scattering_map.q_axis, start_qpoint, end_qpoint)
        return self.advise_experimental_path(
            qpoints,
            coordinates=coordinates,
            preserve_image=preserve_image,
            suggest_supercell=suggest_supercell,
        )

    def plan_map_q_axis(
        self,
        scattering_map,
        start_qpoint,
        end_qpoint,
        coordinates="cartesian",
        q_policy="strict",
        max_error_reduced=None,
        max_error_cartesian=None,
    ):
        """Plan ``Q=q+G`` points from an experimental map's calibrated q axis."""

        if not hasattr(scattering_map, "q_axis"):
            raise TypeError("scattering_map must expose a q_axis attribute")
        qpoints = qpoints_from_path_axis(scattering_map.q_axis, start_qpoint, end_qpoint)
        return self.plan_extended_zone_q_path(
            qpoints,
            coordinates=coordinates,
            q_policy=q_policy,
            max_error_reduced=max_error_reduced,
            max_error_cartesian=max_error_cartesian,
        )

    def build_commensurate_path(self, start_reduced, end_reduced):
        fractions = self._find_allowed_fractions(start_reduced, end_reduced)
        if not fractions:
            return np.zeros((0, 3)), np.zeros((0,), dtype=float)

        start = np.asarray(start_reduced, dtype=float)
        end = np.asarray(end_reduced, dtype=float)
        fvals = np.array([float(f) for f in fractions])
        qpoints = np.array([start + f * (end - start) for f in fvals])
        qpoints = np.where(np.isclose(qpoints, 0.0, atol=1e-14), 0.0, qpoints)
        return qpoints, fvals

    def write_phonopy_qpoints(self, qpoints_reduced, path):
        q = np.atleast_2d(np.asarray(qpoints_reduced, dtype=float))
        if not self.validate_qpoints(q):
            raise ValueError("all phonopy q points must be commensurate with the MD supercell")
        np.savetxt(path, q, fmt="%.10f")

    def _find_allowed_fractions(self, start_reduced, end_reduced):
        start = np.array([self._as_fraction(x) for x in start_reduced])
        end = np.array([self._as_fraction(x) for x in end_reduced])
        if np.all(start == end):
            return [Fraction(0, 1)]

        start_sc = start @ self.P.T
        delta_sc = (end - start) @ self.P.T
        possible = None
        for a, b in zip(start_sc, delta_sc):
            values = self._solve_integer_line(a, b)
            if values is None:
                continue
            if not values:
                return []
            possible = set(values) if possible is None else possible & set(values)
        return sorted(possible) if possible else []

    @staticmethod
    def _solve_integer_line(a, b):
        if b == 0:
            return None if a.denominator == 1 else []

        low, high = (a, a + b) if b > 0 else (a + b, a)
        n_min = math.floor(float(low))
        n_max = math.ceil(float(high))
        vals = [Fraction(n - a, b) for n in range(n_min, n_max + 1)]
        return [x for x in vals if 0 <= x <= 1]

    @staticmethod
    def _as_fraction(value, max_denominator=96, atol=1e-8):
        frac = Fraction(str(float(value)))
        simple = frac.limit_denominator(max_denominator)
        if abs(float(simple) - float(value)) <= atol:
            return simple
        return frac


def suggest_diagonal_supercell(qpoints_reduced, max_denominator=96):
    """Suggest diagonal repeats so every supplied reduced q point is allowed."""

    q = np.atleast_2d(np.asarray(qpoints_reduced, dtype=float))
    repeats = []
    for axis in range(3):
        denom = 1
        for val in q[:, axis]:
            frac = Fraction(str(float(val))).limit_denominator(max_denominator)
            denom = _lcm(denom, frac.denominator)
        repeats.append(denom)
    return np.asarray(repeats, dtype=int)

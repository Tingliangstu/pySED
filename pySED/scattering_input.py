"""Input-file bridge for experiment-facing scattering workflows.

This module keeps ``input_SED.in`` parsing separate from the numerical
scattering APIs.  The first CLI-facing mode is the q-advisor because every
finite-trajectory SED/DSF/EELS calculation needs the same commensurability
decision before spectra are computed.
"""

from pathlib import Path

import numpy as np

from pySED.q_advisor import QAdvisor, qpoints_from_mesh, qpoints_from_path_segments


class ScatteringInputError(ValueError):
    """Error raised for invalid ``input_SED.in`` scattering settings."""

    def __init__(self, message, written_files=None):
        super().__init__(message)
        self.written_files = written_files or {}


def _input_relative_path(params, path_value):
    path = Path(path_value)
    if path.is_absolute():
        return path
    input_path = Path(getattr(params, "input_file", "input_SED.in"))
    return input_path.parent / path


def _output_prefix(params):
    prefix = getattr(params, "scattering_output_prefix", None)
    if prefix:
        return Path(prefix)
    return Path("%s_scattering" % getattr(params, "out_files_name", "pySED"))


def _ensure_output_parent(path):
    parent = Path(path).parent
    if str(parent) not in ("", "."):
        parent.mkdir(parents=True, exist_ok=True)


def build_q_advisor_from_input(params):
    """Build :class:`pySED.q_advisor.QAdvisor` from parsed input parameters."""

    if not hasattr(params, "prim_unitcell"):
        raise ScatteringInputError("scattering mode requires prim_unitcell in input_SED.in")

    primitive = np.asarray(params.prim_unitcell, dtype=float)
    if primitive.shape != (3, 3):
        raise ScatteringInputError("prim_unitcell must contain 9 numeric values")

    supercell_dim = np.asarray(getattr(params, "supercell_dim", [1, 1, 1]), dtype=float)
    if supercell_dim.shape != (3,):
        raise ScatteringInputError("supercell_dim must contain 3 integer repeats")
    if np.any(supercell_dim <= 0):
        raise ScatteringInputError("supercell_dim values must be positive")

    supercell_cell = np.diag(supercell_dim) @ primitive
    return QAdvisor(primitive, supercell_cell)


def build_scattering_qpoints_from_input(params):
    """Return requested Q points and source metadata from parsed input."""

    option = str(getattr(params, "scattering_qpoints_option", "path")).lower()

    if option in ("mesh", "q_mesh", "selected_mesh"):
        mesh = qpoints_from_mesh(
            getattr(params, "scattering_q_mesh_H", 0.0),
            getattr(params, "scattering_q_mesh_K", 0.0),
            getattr(params, "scattering_q_mesh_L", 0.0),
        )
        return mesh.qpoints_reduced, {
            "source": "mesh",
            "mesh_shape": mesh.mesh_shape,
        }

    if option in ("file", "q_file"):
        q_file = getattr(params, "scattering_q_file", None)
        if not q_file:
            raise ScatteringInputError(
                "scattering_qpoints_option = file requires scattering_q_file"
            )
        path = _input_relative_path(params, q_file)
        qpoints = np.loadtxt(path, dtype=float)
        qpoints = np.asarray(qpoints, dtype=float).reshape(-1, 3)
        return qpoints, {"source": "file", "path": str(path)}

    if option in ("points", "inline", "direct"):
        qpoints = getattr(params, "scattering_qpoints", None)
        if qpoints is None:
            raise ScatteringInputError(
                "scattering_qpoints_option = points requires scattering_qpoints"
            )
        return np.asarray(qpoints, dtype=float).reshape(-1, 3), {"source": "points"}

    if option in ("path", "q_path"):
        if not hasattr(params, "q_path"):
            raise ScatteringInputError(
                "scattering_qpoints_option = path requires num_qpaths and q_path"
            )

        vertices = np.asarray(params.q_path, dtype=float)
        if vertices.ndim != 2 or vertices.shape[1] != 3 or vertices.shape[0] < 2:
            raise ScatteringInputError("q_path must contain at least two 3D vertices")

        steps = getattr(params, "scattering_q_path_steps", None)
        if steps is None:
            return vertices, {"source": "path_vertices", "num_segments": vertices.shape[0] - 1}

        steps = np.asarray(steps, dtype=int).ravel()
        num_segments = vertices.shape[0] - 1
        if steps.size == 1:
            steps = np.repeat(steps[0], num_segments)
        if steps.size != num_segments:
            raise ScatteringInputError(
                "scattering_q_path_steps must have one value or num_qpaths values"
            )

        path = qpoints_from_path_segments(
            starts=vertices[:-1],
            ends=vertices[1:],
            steps=steps,
            endpoint=False,
        )
        qpoints = np.vstack([path.qpoints_reduced, vertices[-1].reshape(1, 3)])
        return qpoints, {
            "source": "path",
            "num_segments": num_segments,
            "segment_steps": steps.tolist(),
        }

    raise ScatteringInputError(
        "scattering_qpoints_option must be path, mesh, file, or points"
    )


def _trajectory_frame_count(params):
    total_steps = int(getattr(params, "total_num_steps", 0) or 0)
    stride = int(getattr(params, "output_data_stride", 0) or 0)
    if total_steps > 0 and stride > 0:
        return max(1, total_steps // stride)
    return None


def _enforce_q_policy(report, params):
    policy = str(getattr(params, "scattering_q_policy", "strict")).lower()
    if policy not in ("strict", "nearest"):
        raise ScatteringInputError("scattering_q_policy must be strict or nearest")

    if policy == "strict" and not report.all_commensurate:
        raise ScatteringInputError(
            "non-commensurate Q points found: %d/%d are commensurate; "
            "set scattering_q_policy = nearest to use nearest points, or enlarge "
            "supercell_dim according to the q-advisor summary"
            % (report.num_commensurate, report.num_points)
        )

    max_error_reduced = getattr(params, "scattering_max_error_reduced", None)
    if max_error_reduced is not None and report.max_error_reduced > max_error_reduced:
        raise ScatteringInputError(
            "nearest commensurate Q exceeds scattering_max_error_reduced"
        )

    max_error_cartesian = getattr(params, "scattering_max_error_cartesian", None)
    if max_error_cartesian is not None and report.max_error_cartesian > max_error_cartesian:
        raise ScatteringInputError(
            "nearest commensurate Q exceeds scattering_max_error_cartesian"
        )


def _write_q_advisor_summary(path, report, q_source, efficiency, params):
    lines = [
        "pySED scattering q-advisor summary",
        "",
        "scattering_mode = %s" % getattr(params, "scattering_mode", "q_advisor"),
        "q_source = %s" % q_source.get("source", "unknown"),
        "coordinates = %s" % getattr(params, "scattering_qpoint_coordinates", "reduced"),
        "q_policy = %s" % getattr(params, "scattering_q_policy", "strict"),
        "preserve_image = %s" % bool(getattr(params, "scattering_preserve_image", True)),
        "num_points = %d" % report.num_points,
        "num_commensurate = %d" % report.num_commensurate,
        "all_commensurate = %s" % report.all_commensurate,
        "max_error_reduced = %.12g" % report.max_error_reduced,
        "max_error_cartesian = %.12g" % report.max_error_cartesian,
        "suggested_diagonal_supercell = %s" % report.suggested_diagonal_supercell.tolist(),
        "",
        "Q-point report",
        report.to_table(),
    ]
    if efficiency is not None:
        lines.extend(["", "Selected-Q efficiency", efficiency.to_table()])

    _ensure_output_parent(path)
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_q_advisor_from_input(params):
    """Run q-advisor mode from ``input_SED.in`` and write report artifacts."""

    advisor = build_q_advisor_from_input(params)
    qpoints, q_source = build_scattering_qpoints_from_input(params)
    coordinates = getattr(params, "scattering_qpoint_coordinates", "reduced")
    preserve_image = bool(getattr(params, "scattering_preserve_image", True))
    report = advisor.advise_experimental_path(
        qpoints,
        coordinates=coordinates,
        preserve_image=preserve_image,
    )

    prefix = _output_prefix(params)
    written = {}

    csv_path = prefix.with_name(prefix.name + "_q_advisor.csv")
    _ensure_output_parent(csv_path)
    report.write_csv(csv_path)
    written["q_advisor_csv"] = str(csv_path)

    nearest_path = prefix.with_name(prefix.name + "_nearest_qpoints_reduced.dat")
    _ensure_output_parent(nearest_path)
    np.savetxt(nearest_path, report.nearest_reduced, fmt="%.10f")
    written["nearest_qpoints_reduced"] = str(nearest_path)

    requested_path = prefix.with_name(prefix.name + "_requested_qpoints_reduced.dat")
    _ensure_output_parent(requested_path)
    np.savetxt(requested_path, report.requested_reduced, fmt="%.10f")
    written["requested_qpoints_reduced"] = str(requested_path)

    efficiency = None
    if bool(getattr(params, "scattering_report_efficiency", True)):
        efficiency = advisor.estimate_selected_q_efficiency(
            report.nearest_reduced,
            num_frames=_trajectory_frame_count(params),
            num_atoms=int(getattr(params, "num_atoms", 0) or 0) or None,
        )

    summary_path = prefix.with_name(prefix.name + "_q_advisor_summary.txt")
    _write_q_advisor_summary(summary_path, report, q_source, efficiency, params)
    written["summary"] = str(summary_path)

    try:
        _enforce_q_policy(report, params)
    except ScatteringInputError as exc:
        exc.written_files.update(written)
        raise

    plan = advisor.plan_extended_zone_q_path(
        qpoints,
        coordinates=coordinates,
        q_policy=getattr(params, "scattering_q_policy", "strict"),
        max_error_reduced=getattr(params, "scattering_max_error_reduced", None),
        max_error_cartesian=getattr(params, "scattering_max_error_cartesian", None),
    )

    plan_csv_path = prefix.with_name(prefix.name + "_extended_zone_plan.csv")
    _ensure_output_parent(plan_csv_path)
    plan.write_csv(plan_csv_path)
    written["extended_zone_plan"] = str(plan_csv_path)

    phonopy_path = prefix.with_name(prefix.name + "_phonopy_qpoints_reduced.dat")
    _ensure_output_parent(phonopy_path)
    plan.write_phonopy_qpoints(phonopy_path)
    written["phonopy_qpoints_reduced"] = str(phonopy_path)

    return {
        "advisor": advisor,
        "report": report,
        "q_source": q_source,
        "efficiency": efficiency,
        "extended_zone_plan": plan,
        "written_files": written,
    }


def run_scattering_from_input(params):
    """Dispatch input-file scattering modes."""

    mode = str(getattr(params, "scattering_mode", "q_advisor")).lower()
    if mode in ("q_advisor", "qadvisor", "advisor"):
        return run_q_advisor_from_input(params)

    raise ScatteringInputError(
        "input_SED.in currently supports scattering_mode = q_advisor; "
        "use the Python scattering APIs for DSF/EELS until their CLI modes are added"
    )

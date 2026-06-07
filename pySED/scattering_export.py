"""Export helpers for paper-ready scattering workflow bundles."""

from pathlib import Path
import csv
import json

import numpy as np

from pySED.compare_experiment import (
    comparison_summary,
    map_value_records,
    write_map_comparison_report,
)
from pySED.visibility_diagnostics import diagnose_mode_visibility


def _json_safe(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    return value


def _write_records_csv(path, records):
    records = list(records)
    fieldnames = list(records[0].keys()) if records else []
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if fieldnames:
            writer.writeheader()
            writer.writerows(records)


def _q_advice_records(q_advice):
    records = []
    for index, advice in enumerate(q_advice or []):
        records.append(
            {
                "index": int(index),
                "is_commensurate": bool(advice.is_commensurate),
                "error_reduced": float(advice.error_reduced),
                "error_cartesian": float(advice.error_cartesian),
                "requested_reduced_1": float(advice.requested_reduced[0]),
                "requested_reduced_2": float(advice.requested_reduced[1]),
                "requested_reduced_3": float(advice.requested_reduced[2]),
                "nearest_reduced_1": float(advice.nearest_reduced[0]),
                "nearest_reduced_2": float(advice.nearest_reduced[1]),
                "nearest_reduced_3": float(advice.nearest_reduced[2]),
                "requested_cartesian_1": float(advice.requested_cartesian[0]),
                "requested_cartesian_2": float(advice.requested_cartesian[1]),
                "requested_cartesian_3": float(advice.requested_cartesian[2]),
                "nearest_cartesian_1": float(advice.nearest_cartesian[0]),
                "nearest_cartesian_2": float(advice.nearest_cartesian[1]),
                "nearest_cartesian_3": float(advice.nearest_cartesian[2]),
                "integer_supercell_index_1": int(advice.integer_supercell_index[0]),
                "integer_supercell_index_2": int(advice.integer_supercell_index[1]),
                "integer_supercell_index_3": int(advice.integer_supercell_index[2]),
            }
        )
    return records


def _q_error_source(result):
    q_plan = getattr(result, "q_plan", None)
    if q_plan is not None:
        return q_plan.report
    return getattr(result, "q_advice", None)


def workflow_export_summary(result, diagnostics=None):
    """Return a JSON-ready summary for a scattering workflow result."""

    summary = {
        "workflow_type": type(result).__name__,
        "has_q_plan": getattr(result, "q_plan", None) is not None,
        "has_comparison": getattr(result, "comparison", None) is not None,
        "has_visibility": hasattr(result, "visibility"),
        "has_scattering_map": getattr(result, "scattering_map", None) is not None,
    }
    q_plan = getattr(result, "q_plan", None)
    if q_plan is not None:
        summary["q_plan"] = q_plan.report.summary()
    elif getattr(result, "q_advice", None):
        q_advice = result.q_advice
        summary["q_advice"] = {
            "num_points": len(q_advice),
            "num_commensurate": int(sum(bool(item.is_commensurate) for item in q_advice)),
            "max_error_reduced": float(max((item.error_reduced for item in q_advice), default=0.0)),
            "max_error_cartesian": float(max((item.error_cartesian for item in q_advice), default=0.0)),
        }
    if getattr(result, "eigen_sed", None) is not None:
        summary["eigen_sed_metadata"] = dict(result.eigen_sed.metadata)
    if getattr(result, "dsf", None) is not None:
        summary["dsf_metadata"] = dict(result.dsf.metadata)
    if getattr(result, "partial_dsf", None) is not None:
        summary["partial_dsf_metadata"] = dict(result.partial_dsf.metadata)
    if getattr(result, "current", None) is not None:
        summary["current_metadata"] = dict(result.current.metadata)
    if getattr(result, "scattering_map", None) is not None:
        summary["scattering_map_metadata"] = dict(result.scattering_map.metadata)
        summary["scattering_map_shape"] = list(result.scattering_map.intensity.shape)
    if getattr(result, "comparison", None) is not None:
        summary["comparison"] = comparison_summary(result.comparison)
    if diagnostics is not None:
        summary["visibility_reason_counts"] = diagnostics.reason_counts()
    return _json_safe(summary)


def write_scattering_workflow_bundle(
    result,
    output_dir,
    prefix="workflow",
    atom_labels=None,
    mode_labels=None,
    include_all_visibility=True,
    visibility_relative_threshold=None,
    linecut_q_indices=None,
    linecut_energy_indices=None,
    write_scattering_map=True,
):
    """Write a reproducible bundle for a DSF/EELS/INS/IXS workflow result.

    The bundle gathers q-advisor decisions, folded phonopy q points, simulated
    maps, visibility diagnostics, and experimental comparison tables when those
    objects are present on ``result``.
    """

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    prefix = str(prefix)
    written = {}

    diagnostics = None
    if hasattr(result, "visibility"):
        diagnostics = diagnose_mode_visibility(
            result.visibility,
            q_advice=_q_error_source(result),
            atom_labels=atom_labels,
            mode_labels=mode_labels,
        )
        diagnostics_path = output / ("%s_visibility_diagnostics.csv" % prefix)
        diagnostics.write_csv(
            diagnostics_path,
            include_all=include_all_visibility,
            relative_threshold=visibility_relative_threshold,
        )
        written["visibility_diagnostics"] = str(diagnostics_path)

    q_plan = getattr(result, "q_plan", None)
    if q_plan is not None:
        q_plan_path = output / ("%s_q_plan.csv" % prefix)
        q_plan.write_csv(q_plan_path)
        written["q_plan"] = str(q_plan_path)
        phonopy_path = output / ("%s_phonopy_qpoints.dat" % prefix)
        q_plan.write_phonopy_qpoints(phonopy_path)
        written["phonopy_qpoints"] = str(phonopy_path)
    elif getattr(result, "q_advice", None):
        q_advice_path = output / ("%s_q_advice.csv" % prefix)
        _write_records_csv(q_advice_path, _q_advice_records(result.q_advice))
        written["q_advice"] = str(q_advice_path)

    if write_scattering_map and getattr(result, "scattering_map", None) is not None:
        map_path = output / ("%s_scattering_map.csv" % prefix)
        _write_records_csv(
            map_path,
            map_value_records(result.scattering_map, value_name="intensity"),
        )
        written["scattering_map"] = str(map_path)

    if getattr(result, "comparison", None) is not None:
        comparison_paths = write_map_comparison_report(
            result.comparison,
            output,
            prefix="%s_comparison" % prefix,
            linecut_q_indices=linecut_q_indices,
            linecut_energy_indices=linecut_energy_indices,
        )
        written.update({"comparison_%s" % key: value for key, value in comparison_paths.items()})

    summary = workflow_export_summary(result, diagnostics=diagnostics)
    summary_path = output / ("%s_summary.json" % prefix)
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    written["summary"] = str(summary_path)
    return written

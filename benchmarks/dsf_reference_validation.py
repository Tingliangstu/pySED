"""Reference validation entry point for pySED scattering development."""

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from pySED.validation import (
    build_dsf_validation_report,
    compare_external_dsf_reference,
    export_dsf_validation_case,
    write_validation_report,
)


def _print_report_summary(report):
    internal = report["internal"]
    coherent = internal["checks"]["coherent"]
    coherent_sem = internal["checks"]["coherent_sem"]
    print("pySED direct-vs-correlation max_abs =", coherent["max_abs"])
    print("pySED direct-vs-correlation sem_max_abs =", coherent_sem["max_abs"])
    for package in report["external"]["packages"]:
        print("%s installed = %s" % (package["module"], package["installed"]))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tolerance", type=float, default=1e-12)
    parser.add_argument("--json", dest="json_path", help="Write a structured JSON report")
    parser.add_argument("--export-case", help="Export the synthetic validation case as .npz")
    parser.add_argument("--compare-reference", help="Compare pySED with an external DSF .npz spectrum")
    parser.add_argument("--reference-key", default="coherent", help="Spectrum key in --compare-reference")
    parser.add_argument("--reference-frequency-key", help="Frequency-axis key in --compare-reference")
    parser.add_argument(
        "--candidate-component",
        default="coherent",
        help="pySED DSF result component to compare against the external spectrum",
    )
    parser.add_argument(
        "--normalization",
        default="least_squares",
        choices=("none", "max", "least_squares"),
        help="Normalization used for external spectrum comparison",
    )
    parser.add_argument("--rmse-tolerance", type=float, help="Optional normalized-RMSE pass threshold")
    parser.add_argument("--comparison-json", help="Write external comparison metrics as JSON")
    parser.add_argument(
        "--require-external",
        action="store_true",
        help="Fail if neither dynasor nor pynamic-style package is importable",
    )
    args = parser.parse_args(argv)

    report = build_dsf_validation_report(tolerance=args.tolerance)
    _print_report_summary(report)

    if args.json_path:
        write_validation_report(report, args.json_path)
        print("validation report written =", args.json_path)

    if args.export_case:
        export_dsf_validation_case(args.export_case)
        print("validation case exported =", args.export_case)

    if args.compare_reference:
        comparison = compare_external_dsf_reference(
            args.compare_reference,
            reference_key=args.reference_key,
            reference_frequency_key=args.reference_frequency_key,
            component=args.candidate_component,
            normalization=args.normalization,
            rmse_tolerance=args.rmse_tolerance,
        )
        print("external comparison normalized_rmse =", comparison["normalized_rmse"])
        print("external comparison correlation =", comparison["correlation"])
        print("external comparison scale_factor =", comparison["scale_factor"])
        if args.comparison_json:
            write_validation_report(comparison, args.comparison_json)
            print("external comparison report written =", args.comparison_json)
        if not comparison["passed"]:
            raise SystemExit("external DSF reference comparison failed")

    if args.require_external:
        installed = {
            package["module"]
            for package in report["external"]["packages"]
            if package.get("installed") and package.get("importable", True)
        }
        if not ({"dynasor", "psf", "pynamic_structure_factor", "pynamic"} & installed):
            raise SystemExit("external DSF validation package is required but unavailable")

    if not report["passed"]:
        raise SystemExit("internal DSF reference validation failed")


if __name__ == "__main__":
    main()

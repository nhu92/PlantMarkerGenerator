#!/usr/bin/env python3
"""Run NCBI table2asn against generated FASTA and feature-table files."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def discrepancy_fatals(path: Path) -> tuple[list[str], list[str]]:
    """Split discrepancy FATAL lines into organelle-generic and unexpected groups."""
    if not path.exists():
        return [], []
    lines = {
        line.strip()
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith("FATAL:")
    }
    expected_prefixes = (
        "FATAL: MISSING_PROTEIN_ID:",
        "FATAL: NO_LOCUS_TAGS:",
        "FATAL: SOURCE_QUALS: taxname (all present, all unique)",
    )
    expected = sorted(line for line in lines if line.startswith(expected_prefixes))
    unexpected = sorted(line for line in lines if not line.startswith(expected_prefixes))
    return expected, unexpected


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run a local NCBI table2asn validation pass."
    )
    parser.add_argument("--table2asn", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--feature-table", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    output = args.outdir / "local_validation.sqn"
    for suffix in (".val", ".stats", ".dr", ".gbf", ".sqn"):
        candidate = output.with_suffix(suffix)
        if candidate.exists():
            candidate.unlink()

    command = [
        str(args.table2asn),
        "-i",
        str(args.fasta),
        "-f",
        str(args.feature_table),
        "-t",
        str(args.template),
        "-o",
        str(output),
        "-M",
        "n",
        "-Z",
        "-V",
        "vb",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    (args.outdir / "table2asn.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (args.outdir / "table2asn.stderr.txt").write_text(completed.stderr, encoding="utf-8")

    validation_file = output.with_suffix(".val")
    statistics_file = output.with_suffix(".stats")
    discrepancy_file = output.with_suffix(".dr")
    validation_lines = []
    if validation_file.exists():
        validation_lines = [
            line
            for line in validation_file.read_text(
                encoding="utf-8", errors="replace"
            ).splitlines()
            if line.strip()
        ]
    expected_fatals, unexpected_fatals = discrepancy_fatals(discrepancy_file)
    passed = (
        completed.returncode == 0
        and not validation_lines
        and not unexpected_fatals
    )
    report = {
        "return_code": completed.returncode,
        "validator_report_exists": validation_file.exists(),
        "validator_report_bytes": validation_file.stat().st_size if validation_file.exists() else 0,
        "validator_message_count": len(validation_lines),
        "statistics_report_exists": statistics_file.exists(),
        "discrepancy_report_exists": discrepancy_file.exists(),
        "organelle_expected_generic_fatals": expected_fatals,
        "unexpected_discrepancy_fatals": unexpected_fatals,
        "zero_validator_messages": not validation_lines,
        "status": "passed_organelle_precheck" if passed else "manual_review_required",
        "output_sqn": str(output.resolve()) if output.exists() else "",
    }
    (args.outdir / "validator_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())

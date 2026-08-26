#!/usr/bin/env python3
"""Run input audit, annotation cleanup, and output audit as one workflow."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent


def run_step(name: str, command: list[str], log_dir: Path) -> int:
    """Run one pipeline step and preserve its stdout and stderr."""
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    (log_dir / f"{name}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (log_dir / f"{name}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    return completed.returncode


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run the core chloroplast annotation cleanup pipeline."
    )
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--genbank", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--metadata-csv", type=Path)
    parser.add_argument("--id-prefix", default="CP")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    log_dir = args.outdir / "logs"
    input_audit = args.outdir / "01_input_audit"
    cleaned = args.outdir / "02_cleaned"
    output_audit = args.outdir / "03_output_audit"
    for directory in (log_dir, input_audit, cleaned, output_audit):
        directory.mkdir(parents=True, exist_ok=True)

    common_config = ["--config", str(args.config)] if args.config else []
    input_return = run_step(
        "01_input_audit",
        [
            sys.executable,
            str(SCRIPT_DIR / "audit_annotations.py"),
            "--genbank",
            str(args.genbank),
            "--fasta",
            str(args.fasta),
            "--outdir",
            str(input_audit),
            *common_config,
        ],
        log_dir,
    )

    clean_command = [
        sys.executable,
        str(SCRIPT_DIR / "clean_annotations.py"),
        "--genbank",
        str(args.genbank),
        "--fasta",
        str(args.fasta),
        "--outdir",
        str(cleaned),
        "--id-prefix",
        args.id_prefix,
        *common_config,
    ]
    if args.metadata_csv:
        clean_command.extend(["--metadata-csv", str(args.metadata_csv)])
    clean_return = run_step("02_clean", clean_command, log_dir)
    if clean_return != 0:
        report = {
            "input_audit_return_code": input_return,
            "clean_return_code": clean_return,
            "output_audit_return_code": None,
            "status": "failed_during_cleanup",
        }
        (args.outdir / "pipeline_summary.json").write_text(
            json.dumps(report, indent=2), encoding="utf-8"
        )
        print(json.dumps(report, indent=2))
        return 2

    output_return = run_step(
        "03_output_audit",
        [
            sys.executable,
            str(SCRIPT_DIR / "audit_annotations.py"),
            "--genbank",
            str(cleaned / "cleaned.gb"),
            "--fasta",
            str(cleaned / "submission.fsa"),
            "--outdir",
            str(output_audit),
            *common_config,
        ],
        log_dir,
    )
    output_summary = json.loads((output_audit / "summary.json").read_text(encoding="utf-8"))
    report = {
        "input_audit_return_code": input_return,
        "clean_return_code": clean_return,
        "output_audit_return_code": output_return,
        "output_error_count": output_summary["error_count"],
        "output_warning_count": output_summary["warning_count"],
        "status": "passed" if output_return == 0 else "manual_review_required",
        "submission_fasta": str((cleaned / "submission.fsa").resolve()),
        "submission_feature_table": str((cleaned / "submission.tbl").resolve()),
        "cleaned_genbank": str((cleaned / "cleaned.gb").resolve()),
    }
    (args.outdir / "pipeline_summary.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2))
    return 0 if output_return == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())


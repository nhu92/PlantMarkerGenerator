#!/usr/bin/env python3
"""Optionally validate organism names against the NCBI Taxonomy API."""

from __future__ import annotations

import argparse
import csv
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path


API_ROOT = "https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon"


def query_batch(names: list[str], retries: int = 3) -> dict[str, dict]:
    """Query a small batch with bounded retries for transient HTTP failures."""
    encoded = ",".join(urllib.parse.quote(name, safe="") for name in names)
    url = f"{API_ROOT}/{encoded}/dataset_report"
    for attempt in range(retries):
        try:
            request = urllib.request.Request(
                url,
                headers={"User-Agent": "chloroplast-submission-pipeline/1.0"},
            )
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = json.load(response)
            found: dict[str, dict] = {}
            for report in payload.get("reports", []):
                taxonomy = report.get("taxonomy", {})
                for query in report.get("query", []):
                    found[query] = taxonomy
            return found
        except Exception:
            if attempt + 1 == retries:
                raise
            time.sleep(2**attempt)
    return {}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate mapping-table organism names against NCBI Taxonomy."
    )
    parser.add_argument("--metadata-csv", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--apply-current-names",
        action="store_true",
        help="Also write reviewed_metadata.csv using NCBI current names where found.",
    )
    args = parser.parse_args()

    with args.metadata_csv.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "organism" not in rows[0]:
        raise ValueError("Metadata CSV must contain an organism column")

    names = list(dict.fromkeys(row["organism"].strip() for row in rows if row["organism"].strip()))
    found: dict[str, dict] = {}
    for start in range(0, len(names), 15):
        batch = names[start : start + 15]
        found.update(query_batch(batch))

    report_rows: list[dict[str, str]] = []
    for row in rows:
        submitted = row["organism"].strip()
        taxonomy = found.get(submitted, {})
        current = taxonomy.get("current_scientific_name", {}).get("name", "")
        report_rows.append(
            {
                "output_id": row.get("output_id", ""),
                "input_fasta_id": row.get("input_fasta_id", ""),
                "submitted_name": submitted,
                "ncbi_taxid": str(taxonomy.get("tax_id", "")),
                "ncbi_current_name": current,
                "status": (
                    "exact"
                    if current == submitted
                    else "synonym_or_changed"
                    if current
                    else "not_found"
                ),
            }
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    with (args.outdir / "taxonomy_report.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(report_rows[0]))
        writer.writeheader()
        writer.writerows(report_rows)

    if args.apply_current_names:
        reviewed_rows = []
        for source_row, report_row in zip(rows, report_rows):
            reviewed = dict(source_row)
            if report_row["ncbi_current_name"]:
                reviewed["organism"] = report_row["ncbi_current_name"]
            reviewed_rows.append(reviewed)
        with (args.outdir / "reviewed_metadata.csv").open(
            "w", encoding="utf-8-sig", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=list(reviewed_rows[0]))
            writer.writeheader()
            writer.writerows(reviewed_rows)

    counts = {
        status: sum(row["status"] == status for row in report_rows)
        for status in ("exact", "synonym_or_changed", "not_found")
    }
    (args.outdir / "taxonomy_summary.json").write_text(
        json.dumps(counts, indent=2), encoding="utf-8"
    )
    print(json.dumps(counts, indent=2))
    return 0 if counts["not_found"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())


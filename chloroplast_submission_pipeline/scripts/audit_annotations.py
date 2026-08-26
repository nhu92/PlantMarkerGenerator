#!/usr/bin/env python3
"""Audit chloroplast GenBank annotations for common submission problems."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from Bio import SeqIO

from pipeline_lib import (
    coding_state,
    feature_gene,
    first_qualifier,
    load_config,
    match_fasta_to_genbank,
    read_records,
)


def has_translation_exception(feature, amino_acid: str) -> bool:
    """Check whether a CDS contains a requested transl_except amino acid."""
    needle = f"aa:{amino_acid}".lower()
    return any(
        needle in str(value).replace(" ", "").lower()
        for value in feature.qualifiers.get("transl_except", [])
    )


def conceptual_translation(feature, state: dict[str, object]) -> str:
    """Apply supported translation exceptions to the conceptual protein."""
    protein = str(state["protein"])
    if protein.endswith("*"):
        protein = protein[:-1]
    if has_translation_exception(feature, "TERM") and protein:
        protein = protein[:-1]
    if has_translation_exception(feature, "Met") and protein:
        protein = "M" + protein[1:]
    return protein


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit a chloroplast GenBank file and write CSV/JSON reports."
    )
    parser.add_argument("--genbank", type=Path, required=True)
    parser.add_argument("--fasta", type=Path)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    starts = set(config["valid_start_codons"])
    stops = set(config["valid_stop_codons"])
    records = read_records(args.genbank, "genbank")
    sequence_identity = None
    if args.fasta:
        fasta_records = read_records(args.fasta, "fasta")
        pairs = match_fasta_to_genbank(fasta_records, records)
        sequence_identity = all(
            str(fasta.seq).upper() == str(genbank.seq).upper()
            for fasta, genbank in pairs
        )

    args.outdir.mkdir(parents=True, exist_ok=True)
    issues: list[dict[str, str]] = []
    record_rows: list[dict[str, object]] = []
    issue_counts: Counter[str] = Counter()
    feature_counts: Counter[str] = Counter()

    def report(record_id: str, feature, code: str, severity: str, detail: str) -> None:
        issue_counts[code] += 1
        issues.append(
            {
                "record_id": record_id,
                "feature_type": feature.type if feature is not None else "source",
                "gene": feature_gene(feature) if feature is not None else "",
                "location": str(feature.location) if feature is not None else "",
                "code": code,
                "severity": severity,
                "detail": detail,
            }
        )

    for record in records:
        before = len(issues)
        local_counts: Counter[str] = Counter()
        source = next((feature for feature in record.features if feature.type == "source"), None)
        organism = first_qualifier(source, "organism")
        if source is None:
            report(record.id, None, "MISSING_SOURCE", "ERROR", "No source feature is present")
        elif not organism or " " not in organism:
            report(
                record.id,
                source,
                "INVALID_ORGANISM",
                "ERROR",
                f"Organism is absent or not a binomial/infraspecific name: {organism!r}",
            )

        exact_seen: Counter[tuple[str, str, str]] = Counter()
        for feature in record.features:
            feature_counts[feature.type] += 1
            local_counts[feature.type] += 1
            gene = feature_gene(feature)
            exact_seen[(feature.type, str(feature.location), gene)] += 1

            if feature.type == "tRNA" and not (
                int(config["trna_minimum_nt"])
                <= len(feature.location)
                <= int(config["trna_maximum_nt"])
            ):
                report(
                    record.id,
                    feature,
                    "IMPLAUSIBLE_TRNA_LENGTH",
                    "WARNING",
                    f"tRNA length is {len(feature.location)} nt",
                )

            if feature.type != "CDS":
                continue
            state = coding_state(record, feature)
            if not gene:
                report(record.id, feature, "MISSING_GENE", "ERROR", "CDS has no gene qualifier")
            if not first_qualifier(feature, "product"):
                report(record.id, feature, "MISSING_PRODUCT", "ERROR", "CDS has no product")
            if state["modulo_three"]:
                report(
                    record.id,
                    feature,
                    "CDS_NOT_MOD3",
                    "ERROR",
                    f"Translated nucleotide length is not divisible by three",
                )
            if state["internal_stops"]:
                report(
                    record.id,
                    feature,
                    "INTERNAL_STOP",
                    "ERROR",
                    f"Translation has {state['internal_stops']} internal stop codon(s)",
                )
            if state["start"] not in starts and not has_translation_exception(feature, "Met"):
                report(
                    record.id,
                    feature,
                    "BAD_START",
                    "ERROR",
                    f"Unexpected table-11 start codon {state['start']}",
                )
            if state["stop"] not in stops and not has_translation_exception(feature, "TERM"):
                report(
                    record.id,
                    feature,
                    "NO_STOP",
                    "ERROR",
                    f"Terminal codon {state['stop']} is not a stop codon",
                )
            provided = first_qualifier(feature, "translation").replace(" ", "")
            if provided and provided.rstrip("*") != conceptual_translation(feature, state):
                report(
                    record.id,
                    feature,
                    "TRANSLATION_MISMATCH",
                    "ERROR",
                    "Provided translation differs from the conceptual translation",
                )

        for (feature_type, location, gene), count in exact_seen.items():
            if count > 1:
                duplicate = next(
                    feature
                    for feature in record.features
                    if feature.type == feature_type
                    and str(feature.location) == location
                    and feature_gene(feature) == gene
                )
                report(
                    record.id,
                    duplicate,
                    "EXACT_DUPLICATE_FEATURE",
                    "WARNING",
                    f"Identical feature occurs {count} times",
                )

        record_rows.append(
            {
                "record_id": record.id,
                "length": len(record.seq),
                "organism": organism,
                "gene": local_counts["gene"],
                "CDS": local_counts["CDS"],
                "tRNA": local_counts["tRNA"],
                "rRNA": local_counts["rRNA"],
                "issue_count": len(issues) - before,
            }
        )

    with (args.outdir / "issues.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        fields = ["record_id", "feature_type", "gene", "location", "code", "severity", "detail"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(issues)
    with (args.outdir / "records.csv").open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(record_rows[0]))
        writer.writeheader()
        writer.writerows(record_rows)

    summary = {
        "records": len(records),
        "sequence_bp": sum(len(record.seq) for record in records),
        "sequence_identity_to_fasta": sequence_identity,
        "feature_counts": dict(feature_counts.most_common()),
        "issue_counts": dict(issue_counts.most_common()),
        "issue_total": len(issues),
        "records_with_issues": sum(int(row["issue_count"] > 0) for row in record_rows),
    }
    summary["error_count"] = sum(1 for issue in issues if issue["severity"] == "ERROR")
    summary["warning_count"] = sum(1 for issue in issues if issue["severity"] == "WARNING")
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0 if summary["error_count"] == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())

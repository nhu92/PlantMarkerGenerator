#!/usr/bin/env python3
"""Clean chloroplast GenBank annotations and create NCBI portal inputs."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from copy import deepcopy
from pathlib import Path

from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from pipeline_lib import (
    build_gene_length_medians,
    build_gene_product_consensus,
    clean_cds,
    feature_gene,
    feature_score,
    first_qualifier,
    gene_location,
    is_origin_spanning,
    load_config,
    match_fasta_to_genbank,
    overlap_clusters,
    choose_rnas,
    read_records,
    safe_sequence_id,
    source_metadata,
    write_feature_table_feature,
)


def load_metadata_overrides(path: Path | None) -> dict[str, dict[str, str]]:
    """Load optional user-reviewed metadata keyed by the original FASTA ID."""
    if path is None:
        return {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows or "input_fasta_id" not in rows[0]:
        raise ValueError("Metadata CSV must contain an input_fasta_id column")
    return {row["input_fasta_id"]: row for row in rows}


def apply_metadata_override(metadata: dict[str, str], override: dict[str, str] | None) -> None:
    """Apply only non-empty user values so blank cells do not erase source data."""
    if not override:
        return
    for key in ("organism", "isolate", "collection_date", "geo_loc_name", "specimen_voucher"):
        value = override.get(key, "").strip()
        if value:
            metadata[key] = value


def remove_origin_duplicates(
    features: list[SeqFeature], record_length: int
) -> tuple[list[SeqFeature], list[SeqFeature]]:
    """Remove circular-origin alternatives when an ordinary model also exists."""
    by_gene: dict[str, list[SeqFeature]] = defaultdict(list)
    for feature in features:
        by_gene[feature_gene(feature)].append(feature)
    rejected: list[SeqFeature] = []
    for copies in by_gene.values():
        wrapped = [feature for feature in copies if is_origin_spanning(feature, record_length)]
        ordinary = [feature for feature in copies if feature not in wrapped]
        if wrapped and ordinary:
            rejected.extend(wrapped)
    return [feature for feature in features if feature not in rejected], rejected


def write_fasta(
    records: list[SeqRecord], mappings: list[dict[str, str]], output: Path
) -> None:
    """Write FASTA with source modifiers used by the SP-GenBank organelle route."""
    with output.open("w", encoding="ascii", newline="\n") as handle:
        for record, mapping in zip(records, mappings):
            modifiers = [
                f"[organism={mapping['organism']}]",
                "[location=chloroplast]",
                "[topology=circular]",
                "[completeness=complete]",
                "[gcode=11]",
            ]
            if mapping["isolate"]:
                modifiers.append(f"[isolate={mapping['isolate']}]")
            handle.write(f">{record.id} {' '.join(modifiers)}\n")
            sequence = str(record.seq).upper()
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def write_csv(path: Path, rows: list[dict[str, str]], fallback_fields: list[str]) -> None:
    """Write a UTF-8 CSV with a stable header even when there are no rows."""
    fields = list(rows[0]) if rows else fallback_fields
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Clean chloroplast annotations from matching GenBank and FASTA files."
    )
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--genbank", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--config", type=Path)
    parser.add_argument("--metadata-csv", type=Path)
    parser.add_argument("--id-prefix", default="CP")
    args = parser.parse_args()

    config = load_config(args.config)
    overrides = load_metadata_overrides(args.metadata_csv)
    fasta_records = read_records(args.fasta, "fasta")
    genbank_records = read_records(args.genbank, "genbank")
    pairs = match_fasta_to_genbank(fasta_records, genbank_records)
    medians = build_gene_length_medians(genbank_records, config)
    product_consensus = build_gene_product_consensus(genbank_records)
    args.outdir.mkdir(parents=True, exist_ok=True)

    cleaned_records: list[SeqRecord] = []
    mappings: list[dict[str, str]] = []
    decisions: list[dict[str, str]] = []
    summary: Counter[str] = Counter()
    width = max(3, len(str(len(pairs))))

    for number, (fasta_record, source_record) in enumerate(pairs, 1):
        output_id = safe_sequence_id(args.id_prefix, number, width)
        metadata = source_metadata(fasta_record, source_record)
        apply_metadata_override(metadata, overrides.get(fasta_record.id))
        mapping = {
            "output_id": output_id,
            "input_fasta_id": fasta_record.id,
            "input_genbank_id": source_record.id,
            "original_fasta_description": fasta_record.description,
            **metadata,
        }
        mappings.append(mapping)

        by_gene: dict[str, list[SeqFeature]] = defaultdict(list)
        for feature in source_record.features:
            gene = feature_gene(feature)
            if feature.type == "CDS" and gene and "fragment" not in gene.lower():
                by_gene[gene].append(feature)

        retained_cds: list[SeqFeature] = []
        for gene, models in by_gene.items():
            for cluster in overlap_clusters(models):
                winner = max(
                    cluster,
                    key=lambda feature: feature_score(
                        source_record, feature, medians.get(gene), config
                    ),
                )
                for alternative in cluster:
                    if alternative is winner:
                        continue
                    decisions.append(
                        {
                            "output_id": output_id,
                            "gene": gene,
                            "location": str(alternative.location),
                            "decision": "dropped overlapping alternative model",
                            "annotator": first_qualifier(alternative, "annotator"),
                        }
                    )
                    summary["dropped_alternative_models"] += 1

                product_filled = False
                if not first_qualifier(winner, "product") and gene in product_consensus:
                    winner = deepcopy(winner)
                    winner.qualifiers["product"] = [product_consensus[gene]]
                    product_filled = True

                cleaned, actions = clean_cds(
                    source_record, winner, medians.get(gene), config
                )
                if cleaned is None:
                    decisions.append(
                        {
                            "output_id": output_id,
                            "gene": gene,
                            "location": str(winner.location),
                            "decision": "; ".join(actions),
                            "annotator": first_qualifier(winner, "annotator"),
                        }
                    )
                    summary["dropped_broken_cds"] += 1
                    continue
                if product_filled:
                    decisions.append(
                        {
                            "output_id": output_id,
                            "gene": gene,
                            "location": str(cleaned.location),
                            "decision": "filled product from batch consensus",
                            "annotator": first_qualifier(winner, "annotator"),
                        }
                    )
                    summary["filled_product_from_batch_consensus"] += 1
                retained_cds.append(cleaned)
                summary["retained_cds"] += 1
                for action in actions:
                    decisions.append(
                        {
                            "output_id": output_id,
                            "gene": gene,
                            "location": str(cleaned.location),
                            "decision": action,
                            "annotator": first_qualifier(winner, "annotator"),
                        }
                    )
                    summary[action] += 1

        retained_cds, circular_duplicates = remove_origin_duplicates(
            retained_cds, len(source_record.seq)
        )
        for feature in circular_duplicates:
            summary["retained_cds"] -= 1
            summary["dropped_circular_duplicate_cds"] += 1
            decisions.append(
                {
                    "output_id": output_id,
                    "gene": feature_gene(feature),
                    "location": str(feature.location),
                    "decision": "dropped circular-origin duplicate model",
                    "annotator": "",
                }
            )

        retained_rna = choose_rnas(
            [
                feature
                for feature in source_record.features
                if feature.type in {"tRNA", "rRNA"}
            ],
            config,
        )
        summary["retained_rna"] += len(retained_rna)
        children = sorted(
            retained_cds + retained_rna,
            key=lambda feature: (
                int(feature.location.start),
                feature.type,
                feature_gene(feature),
            ),
        )

        output = SeqRecord(
            source_record.seq,
            id=output_id,
            name=output_id,
            description=f"{metadata['organism']} chloroplast, complete genome",
        )
        output.annotations = {
            "molecule_type": "DNA",
            "topology": "circular",
            "organism": metadata["organism"],
            "source": metadata["organism"],
        }
        source_qualifiers = {
            "organism": [metadata["organism"]],
            "organelle": ["plastid:chloroplast"],
        }
        for key in ("isolate", "collection_date", "geo_loc_name", "specimen_voucher"):
            if metadata[key]:
                source_qualifiers[key] = [metadata[key]]
        output.features = [
            SeqFeature(
                FeatureLocation(0, len(output.seq), strand=1),
                type="source",
                qualifiers=source_qualifiers,
            )
        ]
        for child in children:
            gene_qualifiers = {"gene": [feature_gene(child)]}
            if (
                feature_gene(child) in set(config["trans_spliced_genes"])
                and len(child.location.parts) > 1
            ):
                gene_qualifiers["exception"] = ["trans-splicing"]
            output.features.append(
                SeqFeature(
                    gene_location(child, len(output.seq), config),
                    type="gene",
                    qualifiers=gene_qualifiers,
                )
            )
            output.features.append(child)
        cleaned_records.append(output)

    SeqIO.write(cleaned_records, args.outdir / "cleaned.gb", "genbank")
    write_fasta(cleaned_records, mappings, args.outdir / "submission.fsa")
    with (args.outdir / "submission.tbl").open("w", encoding="utf-8", newline="\n") as handle:
        for record in cleaned_records:
            handle.write(f">Feature {record.id}\n")
            for feature in record.features:
                if feature.type != "source":
                    write_feature_table_feature(handle, feature)

    write_csv(
        args.outdir / "sequence_mapping_and_metadata.csv",
        mappings,
        [
            "output_id",
            "input_fasta_id",
            "input_genbank_id",
            "original_fasta_description",
            "organism",
            "isolate",
            "collection_date",
            "geo_loc_name",
            "specimen_voucher",
        ],
    )
    write_csv(
        args.outdir / "annotation_decisions.csv",
        decisions,
        ["output_id", "gene", "location", "decision", "annotator"],
    )
    summary_payload = {
        "input_records": len(pairs),
        "input_bp": sum(len(record.seq) for record, _ in pairs),
        "output_records": len(cleaned_records),
        "output_bp": sum(len(record.seq) for record in cleaned_records),
        "dna_sequences_unchanged": all(
            str(fasta.seq).upper() == str(cleaned.seq).upper()
            for (fasta, _), cleaned in zip(pairs, cleaned_records)
        ),
        "actions": dict(summary.most_common()),
    }
    (args.outdir / "build_summary.json").write_text(
        json.dumps(summary_payload, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary_payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

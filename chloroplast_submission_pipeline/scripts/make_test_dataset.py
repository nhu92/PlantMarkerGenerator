#!/usr/bin/env python3
"""Create a small exact-record test set from a larger FASTA/GenBank pair."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from Bio import SeqIO

from pipeline_lib import match_fasta_to_genbank, read_records


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Select complete records from matching FASTA and GenBank files."
    )
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--genbank", type=Path, required=True)
    parser.add_argument("--indices", type=int, nargs="+", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()

    fasta_records = read_records(args.fasta, "fasta")
    genbank_records = read_records(args.genbank, "genbank")
    pairs = match_fasta_to_genbank(fasta_records, genbank_records)
    selected = []
    for one_based_index in args.indices:
        if not 1 <= one_based_index <= len(pairs):
            raise IndexError(f"Record index {one_based_index} is outside 1..{len(pairs)}")
        selected.append((one_based_index, pairs[one_based_index - 1]))

    args.outdir.mkdir(parents=True, exist_ok=True)
    SeqIO.write(
        [pair[0] for _, pair in selected],
        args.outdir / "test_input.fasta",
        "fasta",
    )
    SeqIO.write(
        [pair[1] for _, pair in selected],
        args.outdir / "test_input.gb",
        "genbank",
    )
    with (args.outdir / "test_manifest.csv").open(
        "w", encoding="utf-8-sig", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["source_index", "fasta_id", "genbank_id", "length"],
        )
        writer.writeheader()
        for source_index, (fasta_record, genbank_record) in selected:
            writer.writerow(
                {
                    "source_index": source_index,
                    "fasta_id": fasta_record.id,
                    "genbank_id": genbank_record.id,
                    "length": len(fasta_record.seq),
                }
            )
    print(f"Wrote {len(selected)} complete records to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

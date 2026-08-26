#!/usr/bin/env python3
"""Smoke tests for the chloroplast submission pipeline."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = ROOT / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from pipeline_lib import DEFAULT_CONFIG, clean_cds  # noqa: E402


class PipelineSmokeTests(unittest.TestCase):
    """Exercise the bundled data and one non-default reading frame."""

    def test_bundled_records_pass_core_pipeline(self) -> None:
        """Run the complete core workflow and verify its key guarantees."""
        with tempfile.TemporaryDirectory(prefix="chloroplast_pipeline_") as directory:
            outdir = Path(directory) / "result"
            command = [
                sys.executable,
                str(SCRIPT_DIR / "run_pipeline.py"),
                "--fasta",
                str(ROOT / "test_data" / "test_input.fasta"),
                "--genbank",
                str(ROOT / "test_data" / "test_input.gb"),
                "--outdir",
                str(outdir),
                "--config",
                str(ROOT / "config.example.json"),
                "--id-prefix",
                "SMOKE",
            ]
            completed = subprocess.run(command, text=True, capture_output=True, check=False)
            self.assertEqual(completed.returncode, 0, completed.stderr or completed.stdout)

            pipeline_summary = json.loads(
                (outdir / "pipeline_summary.json").read_text(encoding="utf-8")
            )
            build_summary = json.loads(
                (outdir / "02_cleaned" / "build_summary.json").read_text(encoding="utf-8")
            )
            output_audit = json.loads(
                (outdir / "03_output_audit" / "summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(pipeline_summary["status"], "passed")
            self.assertTrue(build_summary["dna_sequences_unchanged"])
            self.assertEqual(build_summary["output_records"], 3)
            self.assertEqual(output_audit["error_count"], 0)
            self.assertEqual(output_audit["warning_count"], 0)
            self.assertTrue((outdir / "02_cleaned" / "submission.tbl").is_file())

            with (ROOT / "test_data" / "test_input.fasta").open() as handle:
                source_sequences = {
                    str(record.seq).upper() for record in SeqIO.parse(handle, "fasta")
                }
            with (outdir / "02_cleaned" / "submission.fsa").open() as handle:
                output_sequences = {
                    str(record.seq).upper() for record in SeqIO.parse(handle, "fasta")
                }
            self.assertEqual(source_sequences, output_sequences)

    def test_codon_start_is_preserved(self) -> None:
        """Ensure a valid offset reading frame is not silently reset to frame one."""
        record = SeqRecord(Seq("AATGAAATAA"), id="offset")
        feature = SeqFeature(
            FeatureLocation(0, 10, strand=1),
            type="CDS",
            qualifiers={
                "gene": ["example"],
                "product": ["example protein"],
                "codon_start": ["2"],
            },
        )
        config = dict(DEFAULT_CONFIG)
        config["minimum_cds_aa"] = 1
        cleaned, actions = clean_cds(record, feature, median_length=2, config=config)
        self.assertIsNotNone(cleaned, actions)
        assert cleaned is not None
        self.assertEqual(cleaned.qualifiers["codon_start"], ["2"])
        self.assertEqual(cleaned.qualifiers["translation"], ["MK"])


if __name__ == "__main__":
    unittest.main(verbosity=2)

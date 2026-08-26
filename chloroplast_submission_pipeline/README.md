# Chloroplast GenBank submission-cleanup pipelines

This folder is a portable, conservative workflow for complete circular chloroplast genomes that already have annotations but fail GenBank checks. The starting inputs are a multi-record GenBank file (`.gb`, `.gbk`, or `.gbff`) and the matching multi-record FASTA file. Geneious files are not required.

The workflow does not change any nucleotide sequence. It ranks competing feature models, preserves credible splicing, applies only configured and documented RNA-editing exceptions, removes models that cannot produce a valid table-11 translation, rebuilds a clean GenBank file, and exports the five-column feature table used by the NCBI organelle submission route.

## Pipeline collection

| Pipeline | Purpose | Entry point | Required? |
|---|---|---|---|
| A | Audit the original annotation | `scripts/audit_annotations.py` | Included in C |
| B | Clean annotations and export GenBank/FASTA/TBL | `scripts/clean_annotations.py` | Included in C |
| C | Run A, B, and a second audit | `scripts/run_pipeline.py` | Yes |
| D | Run a local `table2asn` precheck | `scripts/run_ncbi_validator.py` | Optional |
| 9 | Check organism names against NCBI Taxonomy | `scripts/validate_taxonomy.py` | Optional and separate |

Core flow:

```text
Input FASTA + input GenBank
          |
          v
Exact sequence-content lock (SHA-256)
          |
          v
Input audit -> candidate ranking -> conservative repair/rejection
          |
          v
Cleaned GenBank + submission FASTA + five-column feature table
          |
          v
Output audit -> manual-review queue -> optional NCBI precheck
```

## Requirements

- Python 3.10 or newer
- Biopython 1.88
- NCBI `table2asn` only for optional Pipeline D
- Internet access only for optional Pipeline 9

Create an environment and install the one Python dependency:

```powershell
# Create an isolated environment.
py -3 -m venv .venv

# Activate it in PowerShell.
.\.venv\Scripts\Activate.ps1

# Install the pinned dependency.
python -m pip install -r requirements.txt
```

## Quick test with the bundled records

The `test_data` folder contains three complete records selected from the supplied source collection. They cover ordinary genes, overlapping alternative annotations, an origin-spanning circular model, short/broken CDS calls, and trans-spliced `rps12`. The manifest records their original indices and identifiers.

From this folder, run:

```powershell
# Run the core pipeline on the bundled three-record test set.
python scripts\run_pipeline.py `
  --fasta test_data\test_input.fasta `
  --genbank test_data\test_input.gb `
  --outdir demo_output `
  --config config.example.json `
  --id-prefix DEMO
```

The command succeeds only when the cleaned output audit contains zero errors. The first audit is allowed to return errors because identifying those input defects is its purpose.

Run the automated smoke tests:

```powershell
# Test the full workflow and codon_start preservation.
python tests\test_smoke.py
```

To create another complete-record subset from a larger matching pair:

```powershell
# Select one-based record indices without truncating any sequence.
python scripts\make_test_dataset.py `
  --fasta "D:\data\chloroplasts.fasta" `
  --genbank "D:\data\chloroplasts.gb" `
  --indices 2 11 85 `
  --outdir another_test_set
```

## Run on a real collection

```powershell
# Replace the paths with the matching FASTA and GenBank files.
python scripts\run_pipeline.py `
  --fasta "D:\data\chloroplasts.fasta" `
  --genbank "D:\data\chloroplasts.gb" `
  --outdir "D:\data\chloroplast_submission" `
  --config config.example.json `
  --id-prefix CHL
```

FASTA and GenBank identifiers may differ, but every sequence must have one exact sequence-content match. The pipeline stops instead of guessing if a sequence is missing, duplicated ambiguously, reverse-complemented, rotated, or otherwise different.

## Core decision rules

The defaults target complete land-plant plastomes translated with genetic code 11.

1. Match every FASTA and GenBank record by exact sequence content.
2. Audit missing source metadata, malformed CDS translations, bad starts/stops, duplicate features, and implausible tRNA lengths.
3. For overlapping CDS alternatives, rank intact frame, product/translation evidence, start and stop validity, annotator priority, and agreement with the batch median protein length.
4. Preserve compound CDS locations and mark multi-part `rps12` as trans-spliced.
5. Allow configured ACG-to-AUG RNA-edited starts (`ndhD`, `psbL`, and `rpl2` by default).
6. Allow configured RNA-created terminal stops only when batch length evidence supports them (`petD` by default).
7. Extend a CDS to a nearby in-frame stop only within the configured limit and when length evidence supports the extension.
8. Keep credible tRNA/rRNA models and distinct inverted-repeat copies while removing overlapping alternatives and circular-origin duplicate models.
9. Rebuild source, gene, CDS, tRNA, and rRNA features with submission-safe qualifiers; do not export `protein_id` or `locus_tag` for the SP-GenBank organelle route.
10. Audit the rebuilt file and verify that all DNA sequences are unchanged.

All thresholds and exception whitelists are in `config.example.json`. Review them before applying the pipeline to algae, non-standard genetic codes, linear plastomes, incomplete assemblies, or unusually rearranged chloroplast genomes.

## Outputs

`OUTDIR/02_cleaned` contains the files intended for review or submission:

- `cleaned.gb`: cleaned annotations in GenBank format for Geneious or other viewers.
- `submission.fsa`: unchanged sequences with short safe IDs and source modifiers.
- `submission.tbl`: NCBI five-column feature table for the SP-GenBank organelle Features page.
- `sequence_mapping_and_metadata.csv`: mapping between new IDs and original records plus source metadata.
- `annotation_decisions.csv`: every alternative removal, boundary repair, RNA-editing exception, and rejected CDS.
- `build_summary.json`: record, base-pair, sequence-identity, and action counts.

`OUTDIR/01_input_audit` and `OUTDIR/03_output_audit` contain per-record and per-feature reports. `pipeline_summary.json` is the final machine-readable status.

## Manual correction in Geneious Prime

The pipeline deliberately does not invent a CDS when no defensible model exists. Filter `annotation_decisions.csv` for decisions beginning with `dropped:` and review those cases:

1. Open the original sequence and `cleaned.gb` side by side.
2. Align the questionable gene and translated protein to one or more closely related, curated chloroplast records.
3. Confirm strand, exon order, splice boundaries, and whether the feature crosses coordinate 1.
4. Move the CDS boundary only when an in-frame start/stop and homologous alignment support it.
5. Distinguish a genuine pseudogene or IR-boundary fragment from a damaged full CDS.
6. Export the corrected annotation as GenBank, rerun Pipeline C, and confirm zero output-audit errors.

Geneious is useful for this visual review, but the pipeline itself remains reproducible and does not depend on Geneious.

## Pipeline 9: optional taxonomy validation

Taxonomy is intentionally separated from annotation cleanup because a name change is a biological metadata decision, not a CDS-coordinate repair.

First run Pipeline C. Then query the organisms in its mapping table:

```powershell
# Query NCBI Taxonomy without changing any metadata.
python scripts\validate_taxonomy.py `
  --metadata-csv demo_output\02_cleaned\sequence_mapping_and_metadata.csv `
  --outdir demo_output\09_optional_taxonomy
```

Review `taxonomy_report.csv`. If the current NCBI names and synonym decisions are appropriate, create a proposed override table:

```powershell
# Write proposed current names; review the CSV before reuse.
python scripts\validate_taxonomy.py `
  --metadata-csv demo_output\02_cleaned\sequence_mapping_and_metadata.csv `
  --outdir demo_output\09_optional_taxonomy `
  --apply-current-names
```

After manual review, rerun the core workflow with the approved overrides:

```powershell
# Apply only the reviewed non-empty metadata values.
python scripts\run_pipeline.py `
  --fasta test_data\test_input.fasta `
  --genbank test_data\test_input.gb `
  --metadata-csv demo_output\09_optional_taxonomy\reviewed_metadata.csv `
  --outdir demo_output_taxonomy_reviewed `
  --config config.example.json `
  --id-prefix DEMO
```

## Pipeline D: optional local NCBI precheck

Download the current `table2asn` executable from NCBI and make a real `.sbt` template containing the submitter and publication information. Then run:

```powershell
# Run table2asn with validator and discrepancy-report output.
python scripts\run_ncbi_validator.py `
  --table2asn "D:\tools\table2asn.exe" `
  --fasta demo_output\02_cleaned\submission.fsa `
  --feature-table demo_output\02_cleaned\submission.tbl `
  --template "D:\submission\real_submitter.sbt" `
  --outdir demo_output\04_ncbi_validation
```

For this organelle submission route, NCBI explicitly says to remove `protein_id` and `locus_tag` entries because the portal assigns protein identifiers. Generic `table2asn -Z` discrepancy output may nevertheless label missing protein IDs or locus tags as fatal genome warnings; those generic warnings do not override the organelle-specific instructions. Actual validator messages, unexpected discrepancy categories, and all portal errors must still be resolved.

The generated `.sqn` is suitable only when the `.sbt` contains the real submitter information. The definitive check is uploading `submission.fsa` and `submission.tbl` through the current SP-GenBank organelle workflow and resolving every portal error before submission.

## Scope and limitations

- The workflow repairs annotations, not assemblies. It does not use raw reads or validate structural assembly accuracy.
- An automated zero-error report means the exported syntax and retained feature models passed these checks; it does not guarantee that a curator will agree with every biological boundary.
- Rejected genes remain visible in `annotation_decisions.csv`. Review them rather than assuming that omission is biologically correct.
- Mixed taxonomic groups can have different RNA-editing patterns, introns, genetic codes, and expected gene lengths. Split biologically distant samples into appropriate batches or adjust the configuration.

## Current NCBI references

- [Submitting Mitochondrial and Chloroplast Genomes to GenBank](https://www.ncbi.nlm.nih.gov/genbank/organelle_submit/)
- [Submission of annotation using a five-column feature table](https://www.ncbi.nlm.nih.gov/genbank/feature_table/)
- [table2asn documentation](https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
- [NCBI discrepancy reports](https://www.ncbi.nlm.nih.gov/genbank/asndisc/)

This folder is covered by the repository's MIT license.

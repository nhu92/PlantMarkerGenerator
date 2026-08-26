# Bundled test data

`test_input.fasta` and `test_input.gb` contain the same three complete chloroplast sequences. They were selected without truncation from source-record indices 2, 11, and 85. `test_manifest.csv` records the original identifiers and lengths.

Expected results with `config.example.json`:

- 3 input and output records
- 434,228 input and output base pairs
- DNA sequences unchanged
- output audit: 0 errors and 0 warnings
- 242 retained CDS features
- 124 retained RNA features
- 17 overlapping alternative CDS models rejected
- 7 broken or implausibly short CDS models sent to the manual-review log
- 4 configured RNA-edited starts represented with translation exceptions
- 2 trans-spliced features marked
- 1 circular-origin duplicate CDS model rejected

The local NCBI `table2asn` precheck produced zero validator messages and no unexpected discrepancy FATAL categories. Its generic missing-`protein_id` and missing-`locus_tag` categories are expected because the NCBI organelle instructions explicitly require those qualifiers to be excluded from the uploaded feature table.

Verify the files after copying:

```powershell
# Compare these values with SHA256SUMS.txt.
Get-FileHash test_input.fasta -Algorithm SHA256
Get-FileHash test_input.gb -Algorithm SHA256
Get-FileHash test_manifest.csv -Algorithm SHA256
```


# ChloroMarker

ChloroMarker is an R Shiny application that converts annotated chloroplast genomes into ranked, assay-ready diagnostic marker candidates. It was rebuilt around the workflow described in the GenomeTrakrCP manuscript and is designed for both local use and deployment to shinyapps.io.

## What the app does

1. Accepts one multi-record GenBank file, multiple `.gb/.gbk/.gbff` files, or an NCBI `.fsa + .tbl` annotation pair. An optional taxonomy-validation CSV can update current species names.
2. Extracts annotated genes and adjacent intergenic spacers, including short gene flanks for primer design.
3. Performs a fast k-mer, positional-diversity, length, and coverage screen before running expensive multiple alignments.
4. Ranks loci by coverage, alignment quality, p-distance, barcode gap, and number of distinguishable species.
5. Evaluates the strongest two-marker combinations.
6. Designs shared IUPAC/degenerate PCR primers with configurable amplicon and degeneracy limits.
7. Predicts amplicon sizes, distance matrices, virtual gel bands, restriction sites, and PCR-RFLP profiles.
8. Recommends one of four interpretation paths:

   - PCR plus gel electrophoresis;
   - PCR-RFLP plus gel electrophoresis;
   - PCR plus Sanger/amplicon sequencing;
   - escalation to nuclear markers or orthogonal evidence.

All primer, gel, and restriction results are in-silico predictions and require wet-lab validation.

## Why the implementation is faster

The original application repeatedly read directories, aligned every requested region from scratch, wrote intermediate files, and loaded several unused analysis packages. ChloroMarker 2.0 now:

- parses uploaded data once and keeps it in memory;
- vectorizes NCBI feature-table parsing;
- removes duplicate gene/CDS/product aliases;
- screens all loci cheaply and precisely aligns only the selected candidate set;
- restricts DECIPHER to one processor for shinyapps.io compatibility;
- caches completed analysis settings within the Shiny session;
- performs calculations in memory and writes files only when the user downloads a result;
- removes all example analyses and file-system side effects from sourced production code.

On the supplied 110-record annotation set, `.fsa + .tbl` parsing takes about one second on the development machine. Candidate extraction takes roughly ten seconds once for the complete 110-taxon set and is substantially faster for the intended genus- or family-level subsets.

## Run locally

Install R 4.3 or newer and the required packages:

```r
install.packages(c("shiny", "bslib", "DT", "ggplot2", "digest", "ape", "BiocManager"))
BiocManager::install(c("Biostrings", "DECIPHER"))
```

Then launch from the repository root:

```r
shiny::runApp(".")
```

`Fingerable.R` remains as a backward-compatible launcher, and `chlgene_fingerprintable.R` remains as a backward-compatible function loader.

## Recommended workflow

1. Upload the annotated genomes.
2. Select a biologically meaningful comparison set, usually species within a genus or family. Comparing highly divergent plant lineages can produce alignments that are irrelevant to assay design.
3. Start with locus coverage `0.8`, 16 precisely aligned candidates, and barcode-gap threshold `0.01`.
4. Inspect the ranked single-marker and two-marker tables.
5. Select a marker and design primers.
6. Review the decision recommendation, primer table, predicted gel, and restriction screen.
7. Download the marker table, distance matrix, primers, and RFLP results.

## Prunus validation

The implementation was checked against the five-species case study in the manuscript: *Prunus persica*, *P. caroliniana*, *P. africana*, *P. sibirica*, and *P. padus*. With a 1% barcode-gap threshold, the app ranked `rpl32–trnL(UAG)` as a top marker and distinguished 5/5 species. It generated shared primers around this region and predicted both amplicon-size and restriction-profile differences, reproducing the manuscript's core assay-development logic.

Predicted sizes can differ from a previously chosen primer pair because the app ranks multiple valid primer pairs under the current Tm, GC, degeneracy, and amplicon constraints.

## Output definitions

- **p-distance:** fraction of compared aligned bases that differ, with pairwise deletion of gaps and ambiguous bases.
- **barcode gap:** minimum interspecific distance minus maximum observed intraspecific distance.
- **distinguishable:** barcode gap meets the user-selected threshold.
- **gel-resolved:** predicted length differs from all other taxa under an approximate 2% resolution rule, with a 5 bp minimum.
- **PCR-RFLP resolved:** simulated visible restriction fragments form a unique gel profile after resolution binning.

## Deploy to shinyapps.io

The `.rscignore` file excludes tests, Git metadata, and local scratch files from deployment.

```r
install.packages("rsconnect")
rsconnect::deployApp(appDir = ".", appName = "chloromarker")
```

For large public deployments, keep `max_candidates` conservative and ask users to analyze a genus/family subset. shinyapps.io containers have finite memory and CPU, and multiple sequence alignment remains the most expensive step.

## Repository layout

```text
app.R                         Shiny UI and server
R/io.R                        GenBank, FASTA/FSA, TBL, and taxonomy parsing
R/regions.R                   Gene/intergenic extraction and fast screening
R/analysis.R                  Alignment, distance, ranking, combinations
R/primers.R                   Primer, gel, restriction, recommendation logic
tests/run_tests.R             Parser and core-analysis smoke tests
Fingerable.R                  Legacy launcher
chlgene_fingerprintable.R     Legacy function loader
```

## Scientific limitations

Chloroplast markers may fail in recently diverged groups, hybrids, chloroplast-capture events, introgressed lineages, or comparisons with inadequate intraspecific sampling. A result that is diagnostic for one accession per species is a candidate assay, not proof of species-wide specificity. Add vouchered conspecific samples, close relatives, likely adulterants, and non-target controls before laboratory validation or regulatory use.

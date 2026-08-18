suppressPackageStartupMessages({
  library(Biostrings)
  library(DECIPHER)
  library(ape)
  library(ggplot2)
})

source("R/io.R")
source("R/regions.R")
source("R/analysis.R")
source("R/primers.R")

replace_base <- function(sequence, positions, base) {
  chars <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  chars[positions] <- base
  paste0(chars, collapse = "")
}

format_origin <- function(sequence) {
  starts <- seq(1L, nchar(sequence), by = 60L)
  vapply(starts, function(start) {
    piece <- substr(sequence, start, min(nchar(sequence), start + 59L))
    grouped <- paste(substring(piece, seq(1L, nchar(piece), 10L), pmin(seq(1L, nchar(piece), 10L) + 9L, nchar(piece))), collapse = " ")
    sprintf("%9d %s", start, tolower(grouped))
  }, character(1))
}

make_record <- function(id, organism, sequence) {
  c(
    sprintf("LOCUS       %-12s %d bp DNA circular PLN", id, nchar(sequence)),
    sprintf("DEFINITION  Synthetic %s chloroplast.", organism),
    "FEATURES             Location/Qualifiers",
    "     source          1..1600",
    sprintf("                     /organism=\"%s\"", organism),
    "     gene            101..600",
    "                     /gene=\"rbcL\"",
    "     gene            complement(800..1300)",
    "                     /gene=\"matK\"",
    "ORIGIN",
    format_origin(sequence),
    "//"
  )
}

base_sequence <- paste(rep(c("A", "C", "G", "T"), 400), collapse = "")
variant <- replace_base(base_sequence, seq(121, 561, by = 20), "T")
fixture <- tempfile(fileext = ".gb")
writeLines(c(make_record("TEST01", "Species alpha", base_sequence), make_record("TEST02", "Species beta", variant)), fixture)

records <- parse_genbank_file(fixture)
stopifnot(length(records) == 2L)
stopifnot(all(vapply(records, function(x) nchar(x$sequence), integer(1)) == 1600L))
stopifnot(all(vapply(records, function(x) nrow(representative_gene_features(x$features)), integer(1)) == 2L))

catalog <- build_region_catalog(records, min_coverage = 1)
stopifnot(all(c("gene:rbcL", "gene:matK") %in% names(catalog)))
analysis <- analyze_region(catalog[["gene:rbcL"]], total_records = 2L, threshold = 0.01)
stopifnot(all(dim(analysis$distance) == c(2L, 2L)))
stopifnot(analysis$metrics$mean_distance > 0)

fragments <- digest_amplicon("AAAAGAATTCTTTT", "GAATTC", 1L)
stopifnot(length(fragments) == 2L, sum(fragments) == 14L)
stopifnot(reverse_complement("ACGTR") == "YACGT")

cat("All ChloroMarker tests passed.\n")

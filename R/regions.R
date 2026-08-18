# Candidate gene and intergenic-region extraction.

reverse_complement <- function(sequence) {
  as.character(Biostrings::reverseComplement(Biostrings::DNAString(sequence)))
}

circular_subseq <- function(sequence, start, end) {
  n <- nchar(sequence)
  if (!n || is.na(start) || is.na(end)) return("")
  start <- ((as.integer(start) - 1L) %% n) + 1L
  end <- ((as.integer(end) - 1L) %% n) + 1L
  if (start <= end) substr(sequence, start, end) else paste0(substr(sequence, start, n), substr(sequence, 1L, end))
}

normalize_locus_label <- function(x) {
  x <- trimws(gsub('["\']', "", x))
  x <- gsub("\\s+", "-", x)
  x <- sub("^tRNA-", "trn", x, ignore.case = TRUE)
  x
}

representative_gene_features <- function(features) {
  if (!nrow(features)) return(features)
  features <- features[!is.na(features$label) & nzchar(features$label) & is.finite(features$start) & is.finite(features$end), , drop = FALSE]
  if (!nrow(features)) return(features)
  if (any(tolower(features$type) == "gene")) features <- features[tolower(features$type) == "gene", , drop = FALSE]
  features$label <- vapply(features$label, normalize_locus_label, character(1))
  priority <- match(tolower(features$type), c("gene", "cds", "trna", "rrna"), nomatch = 9L)
  features <- features[order(features$label, priority, -(features$end - features$start)), , drop = FALSE]
  features <- features[!duplicated(features$label), , drop = FALSE]
  features[order(features$start, features$end), , drop = FALSE]
}

extract_record_regions <- function(record, flank_bp = 70L, min_region_length = 40L, max_region_length = 5000L) {
  genome <- record$sequence
  genes <- representative_gene_features(record$features)
  if (!nrow(genes)) return(list())
  out <- list()

  for (i in seq_len(nrow(genes))) {
    feature <- genes[i, ]
    sequence <- circular_subseq(genome, feature$start, feature$end)
    if (feature$strand < 0L) sequence <- reverse_complement(sequence)
    if (nchar(sequence) >= min_region_length && nchar(sequence) <= max_region_length) {
      key <- paste0("gene:", feature$label)
      out[[key]] <- list(
        locus = feature$label, type = "gene", sequence = sequence,
        start = feature$start, end = feature$end, strand = feature$strand
      )
    }
  }

  if (nrow(genes) >= 2L) {
    for (i in seq_len(nrow(genes))) {
      j <- if (i == nrow(genes)) 1L else i + 1L
      left <- genes[i, ]
      right <- genes[j, ]
      genome_length <- nchar(genome)
      gap <- if (j == 1L) (right$start + genome_length) - left$end - 1L else right$start - left$end - 1L
      if (gap < 0L || gap > max_region_length) next
      raw_start <- left$end - flank_bp + 1L
      raw_end <- right$start + flank_bp - 1L
      if (j == 1L) raw_end <- raw_end + genome_length
      sequence <- circular_subseq(genome, raw_start, raw_end)
      labels <- c(left$label, right$label)
      canonical <- sort(labels)
      if (!identical(labels, canonical)) sequence <- reverse_complement(sequence)
      locus <- paste(canonical, collapse = "__")
      if (nchar(sequence) >= min_region_length && nchar(sequence) <= max_region_length + 2L * flank_bp) {
        key <- paste0("igs:", locus)
        candidate <- list(
          locus = locus, type = "intergenic", sequence = sequence,
          start = raw_start, end = raw_end, strand = if (identical(labels, canonical)) 1L else -1L
        )
        if (is.null(out[[key]]) || nchar(candidate$sequence) < nchar(out[[key]]$sequence)) out[[key]] <- candidate
      }
    }
  }
  out
}

build_region_catalog <- function(records, min_coverage = 0.7, flank_bp = 70L, max_region_length = 5000L) {
  catalog <- list()
  for (record in records) {
    regions <- extract_record_regions(record, flank_bp = flank_bp, max_region_length = max_region_length)
    for (key in names(regions)) {
      if (is.null(catalog[[key]])) {
        catalog[[key]] <- list(
          key = key, locus = regions[[key]]$locus, type = regions[[key]]$type,
          sequences = character(), organisms = character()
        )
      }
      catalog[[key]]$sequences[[record$id]] <- regions[[key]]$sequence
      catalog[[key]]$organisms[[record$id]] <- record$organism
    }
  }
  required <- max(2L, ceiling(length(records) * min_coverage))
  catalog[vapply(catalog, function(x) length(x$sequences) >= required, logical(1))]
}

quick_region_metrics <- function(region, total_records) {
  sequences <- unname(region$sequences)
  lengths <- nchar(sequences)
  kmer <- Biostrings::oligonucleotideFrequency(Biostrings::DNAStringSet(sequences), width = 4L, as.prob = TRUE)
  diversity <- if (nrow(kmer) > 1L) mean(stats::dist(kmer)) else 0
  sample_index <- unique(round(seq(1, length(sequences), length.out = min(16L, length(sequences)))))
  sample_sequences <- sequences[sample_index]
  padded_width <- max(nchar(sample_sequences))
  padded <- paste0(sample_sequences, strrep("-", padded_width - nchar(sample_sequences)))
  positional <- if (length(padded) > 1L) {
    mean(Biostrings::stringDist(Biostrings::BStringSet(padded), method = "hamming")) / padded_width
  } else 0
  data.frame(
    key = region$key,
    locus = region$locus,
    type = region$type,
    coverage_n = length(sequences),
    coverage = length(sequences) / total_records,
    median_length = stats::median(lengths),
    min_length = min(lengths),
    max_length = max(lengths),
    length_range = diff(range(lengths)),
    raw_unique = length(unique(sequences)),
    kmer_diversity = diversity,
    positional_diversity = positional,
    quick_score = (diversity + 0.65 * positional) * log1p(stats::median(lengths)) + log1p(diff(range(lengths))) / 10,
    stringsAsFactors = FALSE
  )
}

catalog_summary <- function(catalog, total_records) {
  if (!length(catalog)) return(data.frame())
  out <- do.call(rbind, lapply(catalog, quick_region_metrics, total_records = total_records))
  rownames(out) <- NULL
  out[order(-out$coverage, -out$quick_score, out$median_length), , drop = FALSE]
}

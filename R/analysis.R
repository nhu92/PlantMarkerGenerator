# Alignment, distance, distinguishability, and marker-combination analysis.

align_region <- function(region, processors = 1L) {
  sequences <- Biostrings::DNAStringSet(region$sequences)
  if (length(sequences) < 2L) stop("At least two sequences are required for comparison.")
  aligned <- suppressMessages(DECIPHER::AlignSeqs(
    sequences,
    processors = max(1L, as.integer(processors)),
    verbose = FALSE
  ))
  matrix <- as.matrix(aligned)
  if (is.null(rownames(matrix))) rownames(matrix) <- names(sequences)
  list(strings = aligned, matrix = matrix)
}

alignment_distance <- function(alignment_matrix) {
  dna <- ape::as.DNAbin(alignment_matrix)
  result <- ape::dist.dna(dna, model = "raw", pairwise.deletion = TRUE, as.matrix = TRUE)
  result[!is.finite(result)] <- NA_real_
  diag(result) <- 0
  result
}

species_distance_summary <- function(distance_matrix, organisms, threshold = 0.01) {
  ids <- intersect(rownames(distance_matrix), names(organisms))
  distance_matrix <- distance_matrix[ids, ids, drop = FALSE]
  organisms <- organisms[ids]
  species <- unique(unname(organisms))
  rows <- lapply(species, function(sp) {
    own <- names(organisms)[organisms == sp]
    other <- names(organisms)[organisms != sp]
    within <- if (length(own) >= 2L) distance_matrix[own, own, drop = FALSE][upper.tri(distance_matrix[own, own, drop = FALSE])] else 0
    between <- if (length(other)) as.numeric(distance_matrix[own, other, drop = FALSE]) else NA_real_
    max_within <- if (length(within) && any(is.finite(within))) max(within, na.rm = TRUE) else 0
    min_between <- if (length(between) && any(is.finite(between))) min(between, na.rm = TRUE) else NA_real_
    gap <- min_between - max_within
    data.frame(
      organism = sp,
      accessions = length(own),
      max_intraspecific = max_within,
      min_interspecific = min_between,
      barcode_gap = gap,
      status = if (is.na(min_between)) "insufficient comparison" else if (gap >= threshold) "distinguishable" else "unresolved",
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

alignment_site_metrics <- function(matrix) {
  valid <- matrix %in% c("A", "C", "G", "T", "a", "c", "g", "t")
  gap_fraction <- mean(!valid)
  variable <- apply(matrix, 2L, function(column) {
    bases <- unique(toupper(column[column %in% c("A", "C", "G", "T", "a", "c", "g", "t")]))
    length(bases) > 1L
  })
  list(
    alignment_length = ncol(matrix),
    gap_fraction = gap_fraction,
    variable_sites = sum(variable),
    variable_fraction = mean(variable)
  )
}

analyze_region <- function(region, total_records, threshold = 0.01, processors = 1L) {
  alignment <- align_region(region, processors = processors)
  distance <- alignment_distance(alignment$matrix)
  species <- species_distance_summary(distance, region$organisms, threshold = threshold)
  sites <- alignment_site_metrics(alignment$matrix)
  upper <- distance[upper.tri(distance)]
  upper <- upper[is.finite(upper)]
  metrics <- data.frame(
    key = region$key,
    locus = region$locus,
    type = region$type,
    coverage_n = length(region$sequences),
    coverage = length(region$sequences) / total_records,
    missing_rate = 1 - length(region$sequences) / total_records,
    alignment_length = sites$alignment_length,
    gap_fraction = sites$gap_fraction,
    variable_sites = sites$variable_sites,
    variable_fraction = sites$variable_fraction,
    mean_distance = if (length(upper)) mean(upper) else 0,
    max_distance = if (length(upper)) max(upper) else 0,
    distinguishable_n = sum(species$status == "distinguishable"),
    unresolved_n = sum(species$status == "unresolved"),
    species_n = nrow(species),
    distinguishable_fraction = mean(species$status == "distinguishable"),
    min_barcode_gap = if (any(is.finite(species$barcode_gap))) min(species$barcode_gap, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
  list(region = region, alignment = alignment, distance = distance, species = species, metrics = metrics)
}

candidate_keys <- function(catalog, total_records, max_candidates = 20L) {
  quick <- catalog_summary(catalog, total_records)
  if (!nrow(quick)) return(character())
  quick$priority <- quick$quick_score * (0.5 + 0.5 * quick$coverage)
  quick <- quick[order(-quick$priority, -quick$coverage, quick$median_length), , drop = FALSE]
  traditional <- intersect(c("gene:rbcL", "gene:matK", "igs:psbA__trnH-GUG", "igs:trnL-UAG__trnF-GAA"), quick$key)
  unique(c(traditional, head(quick$key, max_candidates)))
}

analyze_catalog <- function(catalog, total_records, max_candidates = 20L, threshold = 0.01,
                            processors = 1L, progress = NULL) {
  keys <- candidate_keys(catalog, total_records, max_candidates = max_candidates)
  if (!length(keys)) stop("No loci met the requested coverage threshold. Lower minimum coverage or check annotations.")
  analyses <- list()
  errors <- character()
  for (i in seq_along(keys)) {
    key <- keys[[i]]
    if (is.function(progress)) progress(i, length(keys), key)
    result <- tryCatch(
      analyze_region(catalog[[key]], total_records, threshold = threshold, processors = processors),
      error = function(e) e
    )
    if (inherits(result, "error")) errors[[key]] <- conditionMessage(result) else analyses[[key]] <- result
  }
  if (!length(analyses)) stop("All candidate alignments failed. Check sequence and annotation quality.")
  ranking <- do.call(rbind, lapply(analyses, `[[`, "metrics"))
  ranking$rank_score <- ranking$distinguishable_fraction * ranking$coverage * (1 - pmin(ranking$gap_fraction, 0.9))
  ranking <- ranking[order(-ranking$distinguishable_n, -ranking$rank_score, -ranking$mean_distance, ranking$alignment_length), , drop = FALSE]
  rownames(ranking) <- NULL
  analyses <- analyses[ranking$key]
  list(ranking = ranking, analyses = analyses, errors = errors)
}

combine_distance_matrices <- function(analyses) {
  ids <- Reduce(union, lapply(analyses, function(x) rownames(x$distance)))
  sum_matrix <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  count_matrix <- matrix(0L, length(ids), length(ids), dimnames = list(ids, ids))
  for (analysis in analyses) {
    present <- rownames(analysis$distance)
    values <- analysis$distance
    ok <- is.finite(values)
    sum_block <- sum_matrix[present, present, drop = FALSE]
    count_block <- count_matrix[present, present, drop = FALSE]
    sum_block[ok] <- sum_block[ok] + values[ok]
    count_block[ok] <- count_block[ok] + 1L
    sum_matrix[present, present] <- sum_block
    count_matrix[present, present] <- count_block
  }
  out <- sum_matrix / count_matrix
  out[count_matrix == 0L] <- NA_real_
  diag(out) <- 0
  out
}

evaluate_marker_combinations <- function(result, max_k = 2L, pool_size = 8L, threshold = 0.01) {
  keys <- head(result$ranking$key, pool_size)
  if (length(keys) < 2L) return(data.frame())
  organisms <- unlist(unname(lapply(result$analyses[keys], function(x) x$region$organisms)), use.names = TRUE)
  organisms <- organisms[!duplicated(names(organisms))]
  rows <- list()
  for (k in 2L:min(max_k, length(keys))) {
    combinations <- utils::combn(keys, k, simplify = FALSE)
    for (combo in combinations) {
      distance <- combine_distance_matrices(result$analyses[combo])
      species <- species_distance_summary(distance, organisms, threshold = threshold)
      rows[[length(rows) + 1L]] <- data.frame(
        combination = paste(vapply(result$analyses[combo], function(x) x$region$locus, character(1)), collapse = " + "),
        keys = paste(combo, collapse = "|"),
        regions = k,
        distinguishable_n = sum(species$status == "distinguishable"),
        unresolved_n = sum(species$status == "unresolved"),
        species_n = nrow(species),
        distinguishable_fraction = mean(species$status == "distinguishable"),
        mean_distance = mean(distance[upper.tri(distance)], na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  }
  out <- do.call(rbind, rows)
  out[order(-out$distinguishable_n, out$regions, -out$mean_distance), , drop = FALSE]
}

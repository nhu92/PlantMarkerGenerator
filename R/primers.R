# Degenerate primer design, virtual gels, restriction profiles, and assay advice.

iupac_sets <- list(
  A = "A", C = "C", G = "G", T = "T",
  R = c("A", "G"), Y = c("C", "T"), S = c("G", "C"), W = c("A", "T"),
  K = c("G", "T"), M = c("A", "C"), B = c("C", "G", "T"),
  D = c("A", "G", "T"), H = c("A", "C", "T"), V = c("A", "C", "G"),
  N = c("A", "C", "G", "T")
)

iupac_lookup <- setNames(names(iupac_sets), vapply(iupac_sets, function(x) paste(sort(x), collapse = ""), character(1)))

consensus_code <- function(column) {
  bases <- sort(unique(toupper(column[column %in% c("A", "C", "G", "T", "a", "c", "g", "t")])))
  if (!length(bases)) return("N")
  unname(iupac_lookup[paste(bases, collapse = "")]) %||% "N"
}

primer_degeneracy <- function(sequence) {
  codes <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  prod(vapply(codes, function(x) length(iupac_sets[[x]] %||% iupac_sets$N), integer(1)))
}

expected_gc <- function(sequence) {
  codes <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  mean(vapply(codes, function(x) {
    bases <- iupac_sets[[x]] %||% iupac_sets$N
    mean(bases %in% c("G", "C"))
  }, numeric(1)))
}

wallace_tm <- function(sequence) {
  codes <- strsplit(sequence, "", fixed = TRUE)[[1L]]
  sum(vapply(codes, function(x) {
    bases <- iupac_sets[[x]] %||% iupac_sets$N
    mean(ifelse(bases %in% c("G", "C"), 4, 2))
  }, numeric(1)))
}

scan_primer_windows <- function(matrix, lengths = 18:24, max_gap_frequency = 0.1,
                                max_degeneracy = 8, gc_range = c(0.35, 0.65),
                                tm_range = c(50, 68)) {
  width <- ncol(matrix)
  rows <- list()
  for (length_primer in lengths) {
    if (width < length_primer) next
    for (start in seq_len(width - length_primer + 1L)) {
      end <- start + length_primer - 1L
      block <- matrix[, start:end, drop = FALSE]
      gap_frequency <- mean(!(toupper(block) %in% c("A", "C", "G", "T")))
      if (gap_frequency > max_gap_frequency) next
      sequence <- paste0(apply(block, 2L, consensus_code), collapse = "")
      degeneracy <- primer_degeneracy(sequence)
      if (!is.finite(degeneracy) || degeneracy > max_degeneracy) next
      gc <- expected_gc(sequence)
      tm <- wallace_tm(sequence)
      if (gc < gc_range[[1L]] || gc > gc_range[[2L]] || tm < tm_range[[1L]] || tm > tm_range[[2L]]) next
      if (grepl("A{5}|C{5}|G{5}|T{5}", sequence)) next
      three_prime <- substr(sequence, nchar(sequence) - 1L, nchar(sequence))
      rows[[length(rows) + 1L]] <- data.frame(
        start = start, end = end, sequence = sequence, length = length_primer,
        degeneracy = degeneracy, gc = gc, tm = tm, gap_frequency = gap_frequency,
        quality = abs(tm - 58) + 8 * gap_frequency + log2(degeneracy) + ifelse(grepl("[GC]$", three_prime), 0, 0.5),
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, rows)
  out[order(out$quality, out$degeneracy, abs(out$gc - 0.5)), , drop = FALSE]
}

amplicons_for_pair <- function(alignment_matrix, forward_start, reverse_end) {
  block <- alignment_matrix[, forward_start:reverse_end, drop = FALSE]
  sequences <- apply(block, 1L, function(x) paste0(toupper(x[x %in% c("A", "C", "G", "T", "a", "c", "g", "t", "N", "n")]), collapse = ""))
  lengths <- apply(block, 1L, function(x) sum(x != "-"))
  list(sequences = sequences, lengths = lengths)
}

design_primers <- function(analysis, min_amplicon = 120L, max_amplicon = 800L,
                           max_degeneracy = 8L, max_gap_frequency = 0.1, top_n = 10L) {
  matrix <- analysis$alignment$matrix
  candidates <- scan_primer_windows(
    matrix,
    max_gap_frequency = max_gap_frequency,
    max_degeneracy = max_degeneracy
  )
  if (!nrow(candidates)) return(list(table = data.frame(), assays = list()))
  candidates <- head(candidates, 250L)
  pairs <- list()
  for (i in seq_len(nrow(candidates))) {
    forward <- candidates[i, ]
    eligible <- which(
      candidates$start > forward$end &
        candidates$end - forward$start + 1L >= min_amplicon &
        candidates$end - forward$start + 1L <= max_amplicon
    )
    if (!length(eligible)) next
    for (j in head(eligible, 30L)) {
      reverse <- candidates[j, ]
      tm_difference <- abs(forward$tm - reverse$tm)
      if (tm_difference > 5) next
      products <- amplicons_for_pair(matrix, forward$start, reverse$end)
      length_range <- diff(range(products$lengths))
      pair_score <- forward$quality + reverse$quality + tm_difference + abs(stats::median(products$lengths) - 350) / 250 - log1p(length_range) / 3
      pairs[[length(pairs) + 1L]] <- list(
        row = data.frame(
          locus = analysis$region$locus,
          forward_primer = forward$sequence,
          reverse_primer = reverse_complement(reverse$sequence),
          forward_start = forward$start,
          reverse_end = reverse$end,
          forward_tm = round(forward$tm, 1),
          reverse_tm = round(reverse$tm, 1),
          forward_gc = round(100 * forward$gc, 1),
          reverse_gc = round(100 * reverse$gc, 1),
          degeneracy = forward$degeneracy * reverse$degeneracy,
          median_amplicon = stats::median(products$lengths),
          min_amplicon = min(products$lengths),
          max_amplicon = max(products$lengths),
          length_range = length_range,
          score = pair_score,
          stringsAsFactors = FALSE
        ),
        amplicons = products$sequences,
        lengths = products$lengths,
        coordinates = c(forward_start = forward$start, reverse_end = reverse$end)
      )
    }
  }
  if (!length(pairs)) return(list(table = data.frame(), assays = list()))
  scores <- vapply(pairs, function(x) x$row$score, numeric(1))
  pairs <- pairs[order(scores)]
  # Keep distinct primer pairs and favor assays with visible length variation.
  signatures <- vapply(pairs, function(x) paste(x$row$forward_primer, x$row$reverse_primer), character(1))
  pairs <- pairs[!duplicated(signatures)]
  pairs <- head(pairs, top_n)
  table <- do.call(rbind, lapply(seq_along(pairs), function(i) {
    row <- pairs[[i]]$row
    row$assay_id <- paste0("A", i)
    row[, c("assay_id", setdiff(names(row), "assay_id")), drop = FALSE]
  }))
  names(pairs) <- table$assay_id
  list(table = table, assays = pairs)
}

restriction_enzymes <- data.frame(
  enzyme = c("AluI", "BamHI", "DraI", "EcoRI", "HaeIII", "HindIII", "MseI", "MspI", "PstI", "RsaI", "Sau3AI", "SmaI", "TaqI", "XbaI"),
  motif = c("AGCT", "GGATCC", "TTTAAA", "GAATTC", "GGCC", "AAGCTT", "TTAA", "CCGG", "CTGCAG", "GTAC", "GATC", "CCCGGG", "TCGA", "TCTAGA"),
  cut = c(2, 1, 3, 1, 2, 1, 1, 1, 5, 2, 1, 3, 1, 1),
  stringsAsFactors = FALSE
)

digest_amplicon <- function(sequence, motif, cut) {
  sequence <- toupper(gsub("[^ACGT]", "N", sequence))
  hits <- gregexpr(motif, sequence, fixed = TRUE)[[1L]]
  hits <- hits[hits > 0L]
  cuts <- sort(unique(hits + cut - 1L))
  cuts <- cuts[cuts > 0L & cuts < nchar(sequence)]
  sort(diff(c(0L, cuts, nchar(sequence))), decreasing = TRUE)
}

restriction_screen <- function(assay, organisms, min_visible_fragment = 30L) {
  sequences <- assay$amplicons
  ids <- names(sequences)
  rows <- list()
  profiles <- list()
  for (e in seq_len(nrow(restriction_enzymes))) {
    enzyme <- restriction_enzymes[e, ]
    fragments <- lapply(sequences, digest_amplicon, motif = enzyme$motif, cut = enzyme$cut)
    visible <- lapply(fragments, function(x) x[x >= min_visible_fragment])
    signatures <- vapply(visible, function(x) {
      resolution <- ifelse(x < 250, 5, ifelse(x < 500, 10, 20))
      paste(sort(round(x / resolution) * resolution, decreasing = TRUE), collapse = "+")
    }, character(1))
    species_labels <- organisms[ids]
    profile_by_species <- split(signatures, species_labels)
    stable <- vapply(profile_by_species, function(x) length(unique(x)) == 1L, logical(1))
    species_signature <- vapply(profile_by_species, function(x) paste(sort(unique(x)), collapse = "/"), character(1))
    unique_species <- vapply(seq_along(species_signature), function(i) sum(species_signature == species_signature[[i]]) == 1L && stable[[i]], logical(1))
    rows[[e]] <- data.frame(
      enzyme = enzyme$enzyme,
      recognition_site = enzyme$motif,
      distinct_profiles = length(unique(signatures)),
      uniquely_resolved_species = sum(unique_species),
      species_n = length(profile_by_species),
      cut_samples = sum(vapply(fragments, length, integer(1)) > 1L),
      stringsAsFactors = FALSE
    )
    profiles[[enzyme$enzyme]] <- data.frame(
      sequence_id = rep(ids, lengths(visible)),
      organism = rep(unname(organisms[ids]), lengths(visible)),
      fragment_bp = unlist(visible, use.names = FALSE),
      stringsAsFactors = FALSE
    )
  }
  table <- do.call(rbind, rows)
  table <- table[table$distinct_profiles > 1L & table$cut_samples > 0L, , drop = FALSE]
  table <- table[order(-table$uniquely_resolved_species, -table$distinct_profiles, -table$cut_samples), , drop = FALSE]
  list(table = table, profiles = profiles)
}

size_resolution_summary <- function(lengths, organisms) {
  ids <- names(lengths)
  species_lengths <- split(lengths, organisms[ids])
  stable <- vapply(species_lengths, function(x) diff(range(x)) <= max(5, 0.02 * mean(x)), logical(1))
  centers <- vapply(species_lengths, stats::median, numeric(1))
  unique_resolved <- vapply(seq_along(centers), function(i) {
    if (!stable[[i]] || length(centers) < 2L) return(FALSE)
    differences <- abs(centers[[i]] - centers[-i])
    all(differences >= pmax(5, 0.02 * pmin(centers[[i]], centers[-i])))
  }, logical(1))
  list(species_n = length(centers), uniquely_resolved = sum(unique_resolved), centers = centers)
}

assay_recommendation <- function(analysis, primer_result, assay_id = NULL) {
  if (!nrow(primer_result$table)) {
    return(data.frame(method = "Primer design unresolved", detail = "No shared primer pair met the current degeneracy, Tm, and amplicon constraints.", stringsAsFactors = FALSE))
  }
  assay_id <- assay_id %||% primer_result$table$assay_id[[1L]]
  assay <- primer_result$assays[[assay_id]]
  size <- size_resolution_summary(assay$lengths, analysis$region$organisms)
  rflp <- restriction_screen(assay, analysis$region$organisms)
  all_sequence_resolved <- analysis$metrics$distinguishable_n == analysis$metrics$species_n
  if (size$uniquely_resolved == size$species_n) {
    method <- "PCR + gel electrophoresis"
    detail <- "Predicted amplicon sizes uniquely resolve every represented species at an approximate 2% gel-resolution rule."
  } else if (nrow(rflp$table) && rflp$table$uniquely_resolved_species[[1L]] == rflp$table$species_n[[1L]]) {
    method <- "PCR-RFLP + gel electrophoresis"
    detail <- sprintf("Digest the amplicon with %s; predicted restriction profiles uniquely resolve every represented species.", rflp$table$enzyme[[1L]])
  } else if (all_sequence_resolved) {
    method <- "PCR + amplicon sequencing"
    detail <- "Size and tested restriction profiles are not fully diagnostic, but the amplicon sequence is; use Sanger or amplicon sequencing."
  } else {
    method <- "Escalate beyond chloroplast DNA"
    detail <- "This plastid region remains unresolved for at least one species comparison; add nuclear markers or orthogonal evidence."
  }
  data.frame(
    method = method,
    detail = detail,
    gel_unique_species = size$uniquely_resolved,
    species_n = size$species_n,
    best_enzyme = if (nrow(rflp$table)) rflp$table$enzyme[[1L]] else NA_character_,
    stringsAsFactors = FALSE
  )
}

gel_data <- function(assay, organisms, restriction = NULL,
                     ladder = c(50, 100, 150, 200, 300, 400, 500, 700, 1000)) {
  if (is.null(restriction)) {
    bands <- data.frame(
      lane = make.unique(unname(organisms[names(assay$lengths)])),
      fragment_bp = as.numeric(assay$lengths),
      stringsAsFactors = FALSE
    )
  } else {
    screen <- restriction_screen(assay, organisms)
    bands <- screen$profiles[[restriction]]
    lane_ids <- unique(bands$sequence_id)
    lane_map <- stats::setNames(make.unique(unname(organisms[lane_ids])), lane_ids)
    bands$lane <- unname(lane_map[bands$sequence_id])
  }
  marker <- data.frame(lane = "Ladder", fragment_bp = ladder, stringsAsFactors = FALSE)
  rbind(marker, bands[, c("lane", "fragment_bp")])
}

plot_virtual_gel <- function(data) {
  data$lane <- factor(data$lane, levels = unique(data$lane))
  data$y <- -log10(data$fragment_bp)
  ggplot2::ggplot(data, ggplot2::aes(x = lane, y = y, colour = lane == "Ladder")) +
    ggplot2::geom_tile(width = 0.72, height = 0.018, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#45A3FF", `FALSE` = "#F4F7FA")) +
    ggplot2::scale_y_continuous(
      breaks = -log10(sort(unique(data$fragment_bp))),
      labels = sort(unique(data$fragment_bp)),
      name = "Predicted fragment size (bp)"
    ) +
    ggplot2::labs(x = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "#101820", colour = NA),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA)
    )
}

# Input parsing for annotated chloroplast genomes.

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || all(is.na(x))) y else x

clean_dna <- function(x) {
  gsub("[^ACGTRYSWKMBDHVN-]", "", toupper(paste0(x, collapse = "")))
}

header_value <- function(header, key) {
  hit <- regexec(sprintf("\\[%s=([^]]+)\\]", key), header, ignore.case = TRUE)
  value <- regmatches(header, hit)[[1L]]
  if (length(value) >= 2L) value[[2L]] else NA_character_
}

parse_fasta_records <- function(path) {
  dna <- Biostrings::readDNAStringSet(path, format = "fasta", use.names = TRUE)
  headers <- names(dna)
  ids <- sub("\\s.*$", "", headers)
  organisms <- vapply(headers, header_value, character(1), key = "organism")
  organisms[is.na(organisms) | !nzchar(organisms)] <- ids[is.na(organisms) | !nzchar(organisms)]
  sequences <- toupper(as.character(dna))
  stats::setNames(lapply(seq_along(dna), function(i) {
    list(
      id = ids[[i]],
      organism = organisms[[i]],
      sequence = sequences[[i]],
      features = data.frame(),
      source = basename(path)
    )
  }), ids)
}

new_feature <- function(type, location) {
  numbers <- as.integer(regmatches(location, gregexpr("[0-9]+", location))[[1L]])
  list(
    type = type,
    start = if (length(numbers)) min(numbers) else NA_integer_,
    end = if (length(numbers)) max(numbers) else NA_integer_,
    strand = if (grepl("complement", location, fixed = TRUE)) -1L else 1L,
    location = location,
    gene = NA_character_,
    product = NA_character_,
    locus_tag = NA_character_
  )
}

features_to_df <- function(features) {
  if (!length(features)) return(data.frame())
  out <- do.call(rbind, lapply(features, function(x) {
    data.frame(
      type = x$type,
      start = as.integer(x$start),
      end = as.integer(x$end),
      strand = as.integer(x$strand),
      location = x$location,
      gene = x$gene %||% NA_character_,
      product = x$product %||% NA_character_,
      locus_tag = x$locus_tag %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))
  out$label <- ifelse(!is.na(out$gene) & nzchar(out$gene), out$gene,
    ifelse(!is.na(out$locus_tag) & nzchar(out$locus_tag), out$locus_tag, out$product))
  out$label <- trimws(gsub("^tRNA-", "trn", out$label, ignore.case = TRUE))
  out
}

parse_genbank_record <- function(lines, source) {
  locus_line <- grep("^LOCUS\\s+", lines, value = TRUE)[1L]
  id <- if (!is.na(locus_line)) strsplit(trimws(locus_line), "\\s+")[[1L]][2L] else "record"
  origin <- grep("^ORIGIN", lines)[1L]
  sequence <- if (!is.na(origin)) clean_dna(lines[(origin + 1L):length(lines)]) else ""
  definition <- grep("^DEFINITION\\s+", lines, value = TRUE)[1L]
  organism_line <- grep("^  ORGANISM\\s+", lines, value = TRUE)[1L]
  organism <- if (!is.na(organism_line)) sub("^  ORGANISM\\s+", "", organism_line) else NA_character_
  if (is.na(organism) || !nzchar(organism)) {
    source_text <- paste(lines, collapse = " ")
    organism <- sub('.*?/organism="([^"]+)".*', "\\1", source_text)
    if (identical(organism, source_text)) organism <- sub("^DEFINITION\\s+", "", definition %||% id)
  }

  feature_start <- grep("^FEATURES", lines)[1L]
  feature_end <- if (!is.na(origin)) origin - 1L else length(lines)
  features <- list()
  current <- NULL
  current_qualifier <- NULL
  finish <- function() {
    if (!is.null(current)) features[[length(features) + 1L]] <<- current
  }
  if (!is.na(feature_start) && feature_start < feature_end) {
    for (line in lines[(feature_start + 1L):feature_end]) {
      if (grepl("^ {5}\\S", line)) {
        finish()
        match <- regexec("^ {5}(\\S+)\\s+(.+)$", line)
        fields <- regmatches(line, match)[[1L]]
        current <- if (length(fields) >= 3L) new_feature(fields[[2L]], fields[[3L]]) else NULL
        current_qualifier <- NULL
      } else if (!is.null(current) && grepl("^ {21}/", line)) {
        q <- trimws(line)
        qmatch <- regexec('^/([^=]+)=?"?(.*?)"?$', q)
        fields <- regmatches(q, qmatch)[[1L]]
        if (length(fields) >= 3L) {
          current_qualifier <- fields[[2L]]
          value <- sub('"$', "", fields[[3L]])
          if (current_qualifier %in% c("gene", "product", "locus_tag")) current[[current_qualifier]] <- value
        }
      } else if (!is.null(current) && !is.null(current_qualifier) && grepl("^ {21}\\S", line)) {
        if (current_qualifier %in% c("gene", "product", "locus_tag")) {
          current[[current_qualifier]] <- trimws(paste(current[[current_qualifier]], sub('"$', "", trimws(line))))
        }
      }
    }
    finish()
  }
  list(
    id = id,
    organism = trimws(organism),
    sequence = sequence,
    features = features_to_df(features),
    source = source
  )
}

parse_genbank_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  ends <- grep("^//\\s*$", lines)
  if (!length(ends)) ends <- length(lines)
  starts <- c(1L, head(ends, -1L) + 1L)
  records <- lapply(seq_along(ends), function(i) {
    parse_genbank_record(lines[starts[[i]]:ends[[i]]], basename(path))
  })
  ids <- make.unique(vapply(records, `[[`, character(1), "id"))
  for (i in seq_along(records)) records[[i]]$id <- ids[[i]]
  stats::setNames(records, ids)
}

parse_tbl_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  header <- grepl("^>Feature\\s+", lines)
  record_number <- cumsum(header)
  record_ids <- sub("^>Feature\\s+", "", trimws(lines[header]))
  feature_line <- grepl("^[<>]?[0-9]+\\t[<>]?[0-9]+\\t[^\\t]+", lines)
  if (!any(feature_line)) return(stats::setNames(vector("list", length(record_ids)), record_ids))
  feature_number <- cumsum(feature_line)

  feature_parts <- strsplit(lines[feature_line], "\\t", fixed = FALSE)
  feature_df <- data.frame(
    feature_number = feature_number[feature_line],
    record_number = record_number[feature_line],
    first = vapply(feature_parts, function(x) x[[1L]], character(1)),
    second = vapply(feature_parts, function(x) x[[2L]], character(1)),
    type = vapply(feature_parts, function(x) x[[3L]], character(1)),
    stringsAsFactors = FALSE
  )
  # Gene rows already span joined features in NCBI five-column tables. Keeping
  # these rows avoids duplicate CDS/product aliases and is much faster at scale.
  gene_rows <- tolower(feature_df$type) == "gene"
  feature_df <- feature_df[gene_rows, , drop = FALSE]

  qualifier_line <- grepl("^\\t\\t\\t(gene|locus_tag)\\t", lines)
  qualifier_parts <- strsplit(lines[qualifier_line], "\\t", fixed = FALSE)
  qualifier_df <- if (length(qualifier_parts)) data.frame(
    feature_number = feature_number[qualifier_line],
    qualifier = vapply(qualifier_parts, function(x) {
      values <- x[nzchar(x)]
      if (length(values)) values[[1L]] else ""
    }, character(1)),
    value = vapply(qualifier_parts, function(x) {
      values <- x[nzchar(x)]
      if (length(values) >= 2L) paste(values[-1L], collapse = " ") else ""
    }, character(1)),
    stringsAsFactors = FALSE
  ) else data.frame()
  if (nrow(qualifier_df)) {
    qualifier_df <- qualifier_df[qualifier_df$feature_number %in% feature_df$feature_number, , drop = FALSE]
    qualifier_df <- qualifier_df[order(match(qualifier_df$qualifier, c("gene", "locus_tag"))), , drop = FALSE]
    qualifier_df <- qualifier_df[!duplicated(qualifier_df$feature_number), , drop = FALSE]
    labels <- stats::setNames(qualifier_df$value, qualifier_df$feature_number)
    feature_df$gene <- unname(labels[as.character(feature_df$feature_number)])
  } else {
    feature_df$gene <- NA_character_
  }
  feature_df$start <- pmin(as.integer(gsub("[<>]", "", feature_df$first)), as.integer(gsub("[<>]", "", feature_df$second)))
  feature_df$end <- pmax(as.integer(gsub("[<>]", "", feature_df$first)), as.integer(gsub("[<>]", "", feature_df$second)))
  feature_df$strand <- ifelse(as.integer(gsub("[<>]", "", feature_df$first)) > as.integer(gsub("[<>]", "", feature_df$second)), -1L, 1L)
  feature_df$location <- paste(feature_df$first, feature_df$second, sep = "..")
  feature_df$product <- NA_character_
  feature_df$locus_tag <- NA_character_
  feature_df$label <- feature_df$gene

  by_record <- vector("list", length(record_ids))
  names(by_record) <- record_ids
  keep_columns <- c("type", "start", "end", "strand", "location", "gene", "product", "locus_tag", "label")
  split_features <- split(feature_df[, keep_columns, drop = FALSE], feature_df$record_number)
  for (index in names(split_features)) by_record[[record_ids[[as.integer(index)]]]] <- split_features[[index]]
  by_record
}

merge_fasta_tbl <- function(fasta_records, tbl_features) {
  for (id in names(fasta_records)) {
    if (id %in% names(tbl_features)) fasta_records[[id]]$features <- tbl_features[[id]]
  }
  fasta_records
}

read_chloromarker_uploads <- function(upload) {
  stopifnot(is.data.frame(upload), all(c("name", "datapath") %in% names(upload)))
  ext <- tolower(tools::file_ext(upload$name))
  gb_idx <- ext %in% c("gb", "gbk", "gbff", "genbank")
  fasta_idx <- ext %in% c("fa", "fas", "fasta", "fna", "fsa")
  tbl_idx <- ext == "tbl"
  taxonomy_idx <- ext == "csv"
  records <- list()
  if (any(gb_idx)) {
    for (path in upload$datapath[gb_idx]) records <- c(records, parse_genbank_file(path))
  }
  if (any(fasta_idx)) {
    fasta_records <- list()
    for (path in upload$datapath[fasta_idx]) fasta_records <- c(fasta_records, parse_fasta_records(path))
    if (any(tbl_idx)) {
      all_features <- list()
      for (path in upload$datapath[tbl_idx]) all_features <- c(all_features, parse_tbl_file(path))
      fasta_records <- merge_fasta_tbl(fasta_records, all_features)
    }
    records <- c(records, fasta_records)
  }
  if (!length(records)) stop("Upload GenBank files, or an annotated FASTA/FSA file with its feature table (.tbl).")
  ids <- make.unique(vapply(records, `[[`, character(1), "id"))
  for (i in seq_along(records)) records[[i]]$id <- ids[[i]]
  names(records) <- ids
  if (any(taxonomy_idx)) {
    tables <- lapply(upload$datapath[taxonomy_idx], utils::read.csv, stringsAsFactors = FALSE, check.names = FALSE)
    taxonomy <- do.call(rbind, tables)
    if (all(c("sequence_id", "ncbi_current_name") %in% names(taxonomy))) {
      taxonomy <- taxonomy[!duplicated(taxonomy$sequence_id), , drop = FALSE]
      rownames(taxonomy) <- taxonomy$sequence_id
      for (id in intersect(names(records), taxonomy$sequence_id)) {
        name <- taxonomy[id, "ncbi_current_name"]
        if (!is.na(name) && nzchar(name)) records[[id]]$organism <- name
        records[[id]]$taxonomy_status <- if ("status" %in% names(taxonomy)) taxonomy[id, "status"] else NA_character_
        records[[id]]$ncbi_taxid <- if ("ncbi_taxid" %in% names(taxonomy)) as.character(taxonomy[id, "ncbi_taxid"]) else NA_character_
      }
    }
  }
  records
}

dataset_summary <- function(records) {
  data.frame(
    sequence_id = vapply(records, `[[`, character(1), "id"),
    organism = vapply(records, `[[`, character(1), "organism"),
    genome_length = vapply(records, function(x) nchar(x$sequence), integer(1)),
    gc_percent = round(vapply(records, function(x) {
      chars <- strsplit(x$sequence, "", fixed = TRUE)[[1L]]
      100 * sum(chars %in% c("G", "C")) / max(1L, sum(chars %in% c("A", "C", "G", "T")))
    }, numeric(1)), 2),
    annotated_features = vapply(records, function(x) nrow(x$features), integer(1)),
    taxonomy_status = vapply(records, function(x) x$taxonomy_status %||% NA_character_, character(1)),
    ncbi_taxid = vapply(records, function(x) x$ncbi_taxid %||% NA_character_, character(1)),
    source = vapply(records, `[[`, character(1), "source"),
    stringsAsFactors = FALSE
  )
}

# Load required libraries
library(Biostrings)
library(msa)
library(ape)
library(dplyr)
library(tidyr)
library(stringr)
library(combinat)
library(DECIPHER)
library(rprimer)
library(stringdist)
library(phangorn)
library(ggplot2)


# Extract the shortest circular interval between two genes or extract full sequence of one gene
extract_shortest_interval_between_genes <- function(gene1, gene2 = NULL, input_dir = ".", output_fasta) {
  files <- list.files(input_dir, pattern = "\\.gb$", full.names = TRUE)
  fasta_seqs <- list()
  fasta_names <- character()
  
  for (file in files) {
    lines <- readLines(file)
    
    # Extract genome sequence
    origin_start <- grep("^ORIGIN", lines)
    origin_seq <- lines[(origin_start + 1):length(lines)]
    origin_seq <- gsub("[0-9]", "", origin_seq)
    origin_seq <- gsub(" ", "", origin_seq)
    origin_seq <- tolower(paste(origin_seq, collapse = ""))
    genome_len <- nchar(origin_seq)
    
    # Parse features
    feature_lines <- grep("^[[:space:]]{5}(CDS|gene)[[:space:]]", lines)
    gene1_coords <- list()
    gene2_coords <- list()
    
    for (idx in feature_lines) {
      type_line <- lines[idx]
      location_text <- gsub(".*?([0-9,<>()\\.]+)", "\\1", type_line)
      
      qualifiers <- lines[(idx + 1):min(idx + 20, length(lines))]
      gene_val <- sub('.*?/gene="(.*?)".*', "\\1", grep("/gene=", qualifiers, value = TRUE)[1])
      product_val <- sub('.*?/product="(.*?)".*', "\\1", grep("/product=", qualifiers, value = TRUE)[1])
      gene_val[is.na(gene_val)] <- ""
      product_val[is.na(product_val)] <- ""
      
      matches <- gregexpr("[0-9]+\\.\\.[0-9]+", location_text, perl = TRUE)
      joined <- unlist(regmatches(location_text, matches))
      coords_list <- lapply(joined, function(j) as.numeric(strsplit(j, "..", fixed = TRUE)[[1]]))
      
      if (gene_val == gene1 || grepl(gene1, product_val)) {
        gene1_coords <- c(gene1_coords, list(coords_list))
      }
      
      if (!is.null(gene2) && (gene_val == gene2 || grepl(gene2, product_val))) {
        gene2_coords <- c(gene2_coords, coords_list)
      }
    }
    
    if (length(gene1_coords) == 0 || (!is.null(gene2) && length(gene2_coords) == 0)) {
      message(sprintf("⚠️ Missing gene(s) in %s", file))
      next
    }
    
    base_name <- tools::file_path_sans_ext(basename(file))
    
    if (is.null(gene2)) {
      # Merge overlapping or contained blocks
      merged_blocks <- list()
      for (block in gene1_coords) {
        new_range <- range(unlist(block))
        overlap <- FALSE
        for (i in seq_along(merged_blocks)) {
          existing_range <- merged_blocks[[i]]
          if ((new_range[1] >= existing_range[1] && new_range[1] <= existing_range[2]) ||
              (new_range[2] >= existing_range[1] && new_range[2] <= existing_range[2]) ||
              (new_range[1] <= existing_range[1] && new_range[2] >= existing_range[2])) {
            merged_blocks[[i]] <- c(min(existing_range[1], new_range[1]), max(existing_range[2], new_range[2]))
            overlap <- TRUE
            break
          }
        }
        if (!overlap) {
          merged_blocks <- c(merged_blocks, list(new_range))
        }
      }
      
      # Extract and stitch sequences
      seq_blocks <- list()
      for (block in merged_blocks) {
        seq <- substr(origin_seq, block[1], block[2])
        seq_blocks <- c(seq_blocks, list(seq))
      }
      
      final_seq <- paste(seq_blocks, collapse = paste(rep("N", 10), collapse = ""))
      
      if (nchar(final_seq) > 10000) {
        message(sprintf("⚠️ Sequence for %s in %s exceeds 10k bp, skipping.", gene1, file))
        next
      }
      
      header <- sprintf(">%s|gene=%s|copies=%d", base_name, gene1, length(seq_blocks))
      fasta_names <- c(fasta_names, header)
      fasta_seqs <- c(fasta_seqs, final_seq)
    } else {
      # Two-gene mode: find shortest circular interval
      gene1_positions <- sapply(gene1_coords, function(x) mean(x[[1]]))
      gene2_positions <- sapply(gene2_coords, function(x) mean(x))
      
      shortest <- Inf
      best_start <- NULL
      best_end <- NULL
      best_dir <- NULL
      for (p1 in gene1_positions) {
        for (p2 in gene2_positions) {
          d1 <- (p2 - p1 + genome_len) %% genome_len
          d2 <- (p1 - p2 + genome_len) %% genome_len
          if (d1 <= d2 && d1 < shortest) {
            shortest <- d1; best_start <- ceiling(p1); best_end <- ceiling(p2); best_dir <- "forward"
          } else if (d2 < d1 && d2 < shortest) {
            shortest <- d2; best_start <- ceiling(p2); best_end <- ceiling(p1); best_dir <- "reverse"
          }
        }
      }
      
      if (best_start <= best_end) {
        seq <- substr(origin_seq, best_start, best_end)
      } else {
        seq <- paste0(substr(origin_seq, best_start, genome_len), substr(origin_seq, 1, best_end))
      }
      
      if (nchar(seq) > 10000) {
        message(sprintf("⚠️ Sequence from %s to %s in %s exceeds 10k bp, skipping.", gene1, gene2, file))
        next
      }
      
      if (nchar(seq) < 30) {
        message(sprintf("⚠️ Sequence for %s in %s are too short (< 30 bp), skipping.", gene1, gene2, file))
        next
      }
      
      header <- sprintf(">%s|start=%d|end=%d|dir=%s|from=%s|to=%s", base_name, best_start, best_end, best_dir, gene1, gene2)
      fasta_names <- c(fasta_names, header)
      fasta_seqs <- c(fasta_seqs, seq)
    }
  }
  
  if (length(fasta_seqs) > 0) {
    fasta_lines <- unlist(mapply(function(h, s) c(h, s), fasta_names, fasta_seqs, SIMPLIFY = FALSE))
    writeLines(fasta_lines, con = output_fasta)
    message(sprintf("✅ Extracted sequences written to %s", output_fasta))
    return(TRUE)
  } else {
    message("❌ No sequences extracted.")
    return(FALSE)
  }
}

# Function to align fasta sequences and calculate pairwise distances
align_and_calc_distance <- function(input_fasta, output_distances) {
  # Read sequences
  seqs <- readDNAStringSet(input_fasta)
  
  # Multiple sequence alignment
  alignment <- msa(seqs, method = "ClustalW")
  
  # Convert alignment to ape alignment format
  aligned_seqs <- as.DNAbin(alignment)
  
  # Calculate pairwise distance matrix
  dist_matrix <- dist.dna(aligned_seqs, model = "K80")  # Kimura 2-parameter
  
  # Convert matrix to a long-form data frame
  dist_df <- as.data.frame(as.table(as.matrix(dist_matrix)))
  colnames(dist_df) <- c("Sequence1", "Sequence2", "Distance")
  
  # Remove self-comparisons and duplicate comparisons
  dist_df <- dist_df[as.character(dist_df$Sequence1) < as.character(dist_df$Sequence2), ]
  
  # Write distances to a TSV file
  write.table(dist_df, file = output_distances, sep = "\t", row.names = FALSE, quote = FALSE)
  
  message(sprintf("✅ Pairwise distances written to %s", output_distances))
}

classify_species_by_distance <- function(tsv_file, offhit_threshold = 0.6, offhit_ratio = 0.7, indist_threshold = 0.01) {
  df <- read.delim(tsv_file, stringsAsFactors = FALSE)
  
  # Step 1: identify off-hit species
  all_species <- unique(c(df$Sequence1, df$Sequence2))
  offhit_species <- c()
  
  for (sp in all_species) {
    related <- df[df$Sequence1 == sp | df$Sequence2 == sp, ]
    too_far <- sum(related$Distance > offhit_threshold)
    total <- nrow(related)
    if (total > 0 && too_far / total >= offhit_ratio) {
      offhit_species <- c(offhit_species, sp)
    }
  }
  
  # Step 2: filter out off-hit rows
  df_filtered <- df[!(df$Sequence1 %in% offhit_species | df$Sequence2 %in% offhit_species), ]
  kept_species <- unique(c(df_filtered$Sequence1, df_filtered$Sequence2))
  
  # Step 3: identify not-distinguishable species (those involved in any close pair)
  close_pairs <- df_filtered[df_filtered$Distance < indist_threshold, ]
  not_distinguishable <- unique(c(close_pairs$Sequence1, close_pairs$Sequence2))
  
  # Step 4: the rest are distinguishable
  distinguishable <- setdiff(kept_species, not_distinguishable)
  
  # Step 5: print counts and return table
  cat("🧬 Final unique species classification:\n")
  cat("  Off-hit:               ", length(unique(offhit_species)), "\n")
  cat("  Not distinguishable:   ", length(unique(not_distinguishable)), "\n")
  cat("  Distinguishable:       ", length(unique(distinguishable)), "\n")
  
  # Output full classification table
  final_status <- data.frame(
    Species = all_species,
    Status = ifelse(all_species %in% offhit_species, "Off-hit",
                    ifelse(all_species %in% not_distinguishable, "Not distinguishable", "Distinguishable"))
  )
  
  return(final_status)
}

combine_and_evaluate_classification <- function(folder_path, max_k = 3) {
  files <- list.files(folder_path, pattern = "classification.tsv$", full.names = TRUE)
  
  gene_tables <- list()
  useful_genes <- c()
  
  for (file in files) {
    df <- read.delim(file, stringsAsFactors = FALSE)
    df$Species <- sub("\\|.*", "", df$Species)
    gene_name <- gsub("_classification$", "", tools::file_path_sans_ext(basename(file)))
    
    status_counts <- table(df$Status)
    distinguishable_count <- status_counts["Distinguishable"]
    
    if (!is.na(distinguishable_count) && distinguishable_count >= 1) {
      col_to_add <- setNames(df["Status"], gene_name)
      gene_tables[[gene_name]] <- cbind(Species = df$Species, col_to_add)
      useful_genes <- c(useful_genes, gene_name)
    } else {
      message(sprintf("⚠️ Gene %s has no distinguishable power, skipped.", gene_name))
    }
  }
  
  gene_names <- names(gene_tables)
  if (length(gene_names) == 0) {
    message("❌ No useful genes or intergenic regions retained.")
    return(NULL)
  }
  
  # 合并为状态矩阵
  merged_df <- Reduce(function(x, y) merge(x, y, by = "Species", all = TRUE), gene_tables)
  rownames(merged_df) <- merged_df$Species
  merged_df$Species <- NULL
  status_mat <- as.matrix(merged_df)
  N_total <- nrow(status_mat)
  
  # 分类逻辑
  classify_combination <- function(mat) {
    apply(mat, 1, function(statuses) {
      statuses <- statuses[!is.na(statuses)]
      if ("Distinguishable" %in% statuses) {
        return("Distinguishable")
      } else if ("Not distinguishable" %in% statuses) {
        return("Not distinguishable")
      } else {
        return("Off-hit")
      }
    })
  }
  
  results <- data.frame()
  stop_after_next_round <- FALSE
  skip_genes <- c()  # 全鉴别单区块列表
  
  for (k in 1:max_k) {
    message(sprintf("🔍 Testing combinations of size k = %d", k))
    
    if (length(gene_names) < k) break
    combos <- combn(gene_names, k, simplify = FALSE)
    found_full_combo <- FALSE
    
    for (combo in combos) {
      # 跳过包含已全鉴别单区块的组合
      if (k > 1 && any(combo %in% skip_genes)) next
      
      combo_label <- paste(combo, collapse = " + ")
      subset_status <- status_mat[, combo, drop = FALSE]
      
      if (all(is.na(subset_status))) next
      
      final_status <- classify_combination(subset_status)
      counts <- table(factor(final_status, levels = c("Distinguishable", "Not distinguishable", "Off-hit")))
      
      # 计算平均距离（来自组合成员的 distance 文件）
      distances <- c()
      for (g in combo) {
        tag <- gsub(" + ", "__", g)
        dist_file <- file.path(folder_path, paste0(tag, "_distance.tsv"))
        if (file.exists(dist_file)) {
          dist_df <- tryCatch(read.delim(dist_file), error = function(e) NULL)
          if (!is.null(dist_df) && "Distance" %in% colnames(dist_df)) {
            distances <- c(distances, dist_df$Distance)
          }
        }
      }
      avg_dist <- if (length(distances) > 0) mean(distances, na.rm = TRUE) else NA
      
      results <- rbind(results, data.frame(
        Combination = combo_label,
        NumRegions = length(combo),
        Distinguishable = counts["Distinguishable"],
        Not_distinguishable = counts["Not distinguishable"],
        Off_hit = counts["Off-hit"],
        AverageDistance = avg_dist
      ))
      
      if (!is.na(counts["Distinguishable"]) && counts["Distinguishable"] == N_total) {
        found_full_combo <- TRUE
        if (k == 1) skip_genes <- c(skip_genes, combo)  # 标记为“已能单独鉴别”的单区块
      }
    }
    
    if (stop_after_next_round) {
      message(sprintf("✅ Stopping after extra round k = %d", k))
      break
    }
    
    if (found_full_combo) {
      stop_after_next_round <- TRUE
    }
  }
  
  # 排序：先 Distinguishable 多 → 区块数少 → 平均距离大
  results <- results[order(-results$Distinguishable, results$NumRegions, -results$AverageDistance), ]
  
  cat("📊 Summary for selected combinations:\n")
  print(results, row.names = FALSE)
  
  return(results)
}

get_common_genes <- function(input_dir = ".") {
  files <- list.files(input_dir, pattern = "\\.gb$", full.names = TRUE)
  gene_sets <- list()
  
  for (file in files) {
    lines <- readLines(file)
    feature_lines <- grep("^[[:space:]]{5}(CDS|gene)[[:space:]]", lines)
    genes <- c()
    
    for (idx in feature_lines) {
      qualifiers <- lines[(idx + 1):min(idx + 20, length(lines))]
      gene_line <- grep('/gene="', qualifiers, value = TRUE)
      product_line <- grep('/product="', qualifiers, value = TRUE)
      
      gene_val <- if (length(gene_line) > 0) sub('.*?/gene="(.*?)".*', "\\1", gene_line[1]) else ""
      product_val <- if (length(product_line) > 0) sub('.*?/product="(.*?)".*', "\\1", product_line[1]) else ""
      
      if (gene_val != "") {
        genes <- c(genes, gene_val)
      } else if (product_val != "") {
        genes <- c(genes, product_val)
      }
    }
    
    gene_sets[[file]] <- unique(genes)
  }
  
  # intersection
  common_genes <- Reduce(intersect, gene_sets)
  return(common_genes)
}

get_common_intergenic_regions <- function(input_dir = ".") {
  files <- list.files(input_dir, pattern = "\\.gb$", full.names = TRUE)
  region_sets <- list()
  
  for (file in files) {
    lines <- readLines(file)
    feature_lines <- grep("^[[:space:]]{5}(CDS|gene)[[:space:]]", lines)
    
    gene_coords <- data.frame(Gene = character(), Start = numeric(), End = numeric(), stringsAsFactors = FALSE)
    
    for (idx in feature_lines) {
      qualifiers <- lines[(idx + 1):min(idx + 20, length(lines))]
      gene_line <- grep('/gene="', qualifiers, value = TRUE)
      product_line <- grep('/product="', qualifiers, value = TRUE)
      gene_val <- if (length(gene_line) > 0) sub('.*?/gene="(.*?)".*', "\\1", gene_line[1]) else ""
      if (gene_val == "" && length(product_line) > 0) {
        gene_val <- sub('.*?/product="(.*?)".*', "\\1", product_line[1])
      }
      
      location_text <- gsub(".*?([0-9,<>()\\.]+)", "\\1", lines[idx])
      match <- regmatches(location_text, regexpr("[0-9]+\\.\\.[0-9]+", location_text))
      if (length(match) == 0) next
      coords <- as.numeric(unlist(strsplit(match, "\\.\\.")))
      
      if (!is.na(gene_val) && length(coords) == 2) {
        gene_coords <- rbind(gene_coords, data.frame(Gene = gene_val, Start = coords[1], End = coords[2]))
      }
    }
    
    # 按起始位置排序
    gene_coords <- gene_coords[order(gene_coords$Start), ]
    
    # 生成相邻基因对名
    regions <- c()
    for (i in 1:(nrow(gene_coords) - 1)) {
      gene1 <- gene_coords$Gene[i]
      gene2 <- gene_coords$Gene[i + 1]
      regions <- c(regions, paste0(gene1, "__", gene2))
    }
    
    region_sets[[file]] <- unique(regions)
  }
  
  # 求交集
  common_regions <- Reduce(intersect, region_sets)
  return(common_regions)
}

find_visible_indel_windows <- function(fasta_file, window_size = 2500, overlap_size = 10, min_length_diff = 10) {
  library(msa)
  library(Biostrings)
  
  seqs <- readDNAStringSet(fasta_file)
  if (length(seqs) < 2) {
    message(sprintf("⚠️ Not enough sequences in %s for comparison.", fasta_file))
    return(NULL)
  }
  
  alignment <- msa(seqs, method = "ClustalW")
  aln_mat <- as.matrix(alignment)
  
  nseq <- nrow(aln_mat)
  seqlen <- ncol(aln_mat)
  
  if (seqlen < 100) {
    message(sprintf("❌ Aligned sequence too short (%d bp) in %s, skipping.", seqlen, fasta_file))
    return(NULL)
  }
  
  step <- window_size - overlap_size
  starts <- seq(1, max(1, seqlen - window_size + 1), by = step)
  results <- data.frame()
  base_name <- tools::file_path_sans_ext(basename(fasta_file))
  out_dir <- dirname(fasta_file)
  
  win_id <- 1
  for (start in starts) {
    end <- min(start + window_size - 1, seqlen)
    submat <- aln_mat[, start:end, drop = FALSE]
    
    gap_counts <- apply(submat, 1, function(row) sum(row == '-'))
    effective_lengths <- (end - start + 1) - gap_counts
    
    len_diff <- max(effective_lengths) - min(effective_lengths)
    
    if (len_diff >= min_length_diff) {
      # 保存 summary 记录
      results <- rbind(results, data.frame(
        Start = start,
        End = end,
        WindowSize = end - start + 1,
        MaxDiff = len_diff
      ))
      
      # 保存详细窗口内容
      species_names <- rownames(submat)
      seqs_window <- apply(submat, 1, paste0, collapse = "")
      out_df <- data.frame(
        Species = species_names,
        SequenceWindow = seqs_window,
        EffectiveLength = effective_lengths,
        stringsAsFactors = FALSE
      )
      
      out_file <- file.path(out_dir, paste0(base_name, "_visibleindel_detail_window", win_id, ".tsv"))
      write.table(out_df, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
      win_id <- win_id + 1
    }
  }
  
  if (nrow(results) == 0) {
    message(sprintf("❌ No visible indel found in %s", fasta_file))
    return(NULL)
  } else {
    message(sprintf("✅ Found %d visible indel windows in %s", nrow(results), fasta_file))
    return(results)
  }
}

process_all_common_regions <- function(input_dir = ".", output_dir = ".", t_offhit = 0.6, r_offhit = 0.7, t_indist = 0.01, mode = c("distance", "length_indel")) {
  
  mode <- match.arg(mode)
  # 1️⃣ 提取共有基因 & 共有基因间隔对
  common_genes <- get_common_genes(input_dir)
  common_intergenic_regions <- get_common_intergenic_regions(input_dir)

  # 2️⃣ 合并两者作为统一的候选区域列表
  all_targets <- c(common_genes, common_intergenic_regions)

  for (target in all_targets) {
    # 识别是基因还是间隔对
    if (grepl("__", target)) {
      parts <- strsplit(target, "__")[[1]]
      g1 <- parts[1]
      g2 <- parts[2]
      tag <- paste0(g1, "__", g2)
      fasta_file <- file.path(output_dir, paste0(tag, ".fasta"))
      dist_file <- file.path(output_dir, paste0(tag, "_distance.tsv"))
      class_file <- file.path(output_dir, paste0(tag, "_classification.tsv"))

      success <- extract_shortest_interval_between_genes(g1, g2, input_dir = input_dir, output_fasta = fasta_file)
    } else {
      g1 <- target
      tag <- g1
      fasta_file <- file.path(output_dir, paste0(tag, ".fasta"))
      dist_file <- file.path(output_dir, paste0(tag, "_distance.tsv"))
      class_file <- file.path(output_dir, paste0(tag, "_classification.tsv"))

      success <- extract_shortest_interval_between_genes(g1, input_dir = input_dir, output_fasta = fasta_file)
    }

    if (!success) {
      message(sprintf("⚠️ Skipping %s due to extraction failure.\n", target))
      next
    }

    if (mode == "distance") {
      align_and_calc_distance(fasta_file, dist_file)
      result <- classify_species_by_distance(dist_file,
                                             offhit_threshold = t_offhit,
                                             offhit_ratio = r_offhit,
                                             indist_threshold = t_indist)
      write.table(result, file = class_file, sep = "\t", row.names = FALSE, quote = FALSE)
      
    } else if (mode == "length_indel") {
      indel_result <- find_visible_indel_windows(fasta_file)
      if (!is.null(indel_result)) {
        write.table(indel_result, file = file.path(output_dir, paste0(tag, "_visibleindel.tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        # Try PCR primer design and simulate
        screen_fasta_for_pcr_indels(fasta_file, output_dir = output_dir)
      }
    }
  }

  # 3️⃣ 汇总筛选出有鉴别力的区域，组合分析
  if (mode == "distance") {
    summary <- combine_and_evaluate_classification(output_dir)
    write.table(summary, file = file.path(output_dir, "combined_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    cat("🎉 All regions processed and evaluated.\n")
  } else {
    cat("🎉 All regions scanned for visible indels.\n")
  }
}

screen_fasta_for_pcr_indels <- function(
    fasta_file,
    min_diff   = 10,
    output_dir = ".",
    log_file   = "failed_primer_design.txt",
    max_mm_amp = 2,          # max mismatches allowed in AmplifyDNA
    min_psize  = 50,
    max_psize  = 2500,
    n_sets     = 3) {
  
  cat("🔬  Screening:", fasta_file, "\n")
  base <- tools::file_path_sans_ext(basename(fasta_file))
  
  # 1. 读取序列
  seqs <- readDNAStringSet(fasta_file)
  if (length(seqs) < 2) {
    msg <- sprintf("⚠️  <2 sequences in %s\n", fasta_file)
    cat(msg); write(msg, file = log_file, append = TRUE)
    return(invisible(NULL))
  }
  
  # 2. 如有 gap 说明已比对；否则先比对
  have_gap <- any(grepl("-", as.character(seqs)))
  aln <- if (have_gap) {
    # 确保所有序列等长，否则报错
    lens <- unique(width(seqs))
    if (length(lens) != 1)
      stop("Sequences contain gaps but are unequal length – not a valid MSA.")
    seqs
  } else {
    AlignSeqs(seqs)
  }
  
  # 3. 引物设计（带 tryCatch）
  primers <- tryCatch(
    DesignPrimers(aln,
                  minProductSize = min_psize,
                  maxProductSize = max_psize,
                  numPrimerSets  = n_sets),
    error = function(e) {
      msg <- sprintf("❌  Primer design failed for %s: %s\n", fasta_file, e$message)
      cat(msg); write(msg, file = log_file, append = TRUE)
      return(NULL)
    })
  
  if (is.null(primers) || length(primers) == 0) {
    msg <- sprintf("⚠️  No primers returned for %s\n", fasta_file)
    cat(msg); write(msg, file = log_file, append = TRUE)
    return(invisible(NULL))
  }
  
  # 4. 模拟扩增
  amps <- AmplifyDNA(aln, primers, maxMismatches = max_mm_amp)
  if (length(amps) == 0) {
    msg <- sprintf("❌  No amplicons simulated for %s\n", fasta_file)
    cat(msg); write(msg, file = log_file, append = TRUE)
    return(invisible(NULL))
  }
  
  # 5. 长度差判定
  lens <- width(amps); names(lens) <- names(amps)
  diff_mat <- outer(lens, lens, function(a,b) abs(a-b)); diag(diff_mat) <- 0
  if (!any(diff_mat >= min_diff)) {
    msg <- sprintf("⚠️  No amplicon-size difference ≥%dbp in %s\n", min_diff, fasta_file)
    cat(msg); write(msg, file = log_file, append = TRUE)
    return(invisible(NULL))
  }
  
  # 6. 导出结果
  primer_tsv <- file.path(output_dir, paste0(base, "_pcr_primers.tsv"))
  amp_tsv    <- file.path(output_dir, paste0(base, "_pcr_amplicons.tsv"))
  write.table(as.data.frame(primers), primer_tsv, sep="\t", row.names=FALSE, quote=FALSE)
  write.table(data.frame(Species = names(lens), AmpliconLength = lens),
              amp_tsv, sep="\t", row.names=FALSE, quote=FALSE)
  cat("✅  Primers saved to", primer_tsv, "\n")
  invisible(list(primers=primers, lengths=lens))
}

findAndDesignGappyAmplicons <- function(
    fasta_file,
    window_size = 100,        # Window size for gappy search
    window_step = 20,         # Sliding step
    top_n_windows = 3,        # How many gappy regions to try
    min_amplicon = 50, max_amplicon = 600,
    overlap_fraction = 0.7, max_per_amp_region = 3,
    design_args = list(),
    output_top = 20           # Maximum number in final output
) {
  # 1. Read raw sequences
  raw_dna <- readDNAStringSet(fasta_file)
  # 2. Perform multiple sequence alignment (e.g. ClustalW)
  alignment <- msa(raw_dna, method = "ClustalW")
  # 3. Convert to DECIPHER format for downstream use
  aln_msa <- DNAMultipleAlignment(alignment)
  aln_mat <- as.matrix(aln_msa)
  aln_len <- ncol(aln_mat)
  seq_names <- rownames(aln_mat)
  
  # --- 1. Sliding window to find gappy regions ---
  window_info <- data.frame()
  for (start in seq(1, aln_len - window_size + 1, by = window_step)) {
    end <- start + window_size - 1
    win <- aln_mat[, start:end, drop=FALSE]
    gap_prop <- mean(win == "-")
    window_info <- rbind(window_info,
                         data.frame(start = start, end = end, size = window_size, gap_prop = gap_prop)
    )
  }
  gappy_regions <- window_info %>%
    arrange(desc(gap_prop)) %>%
    slice_head(n = top_n_windows)
  
  if (nrow(gappy_regions) == 0) {
    message("No gappy regions found in alignment.")
    return(NULL)
  }
  
  # --- 2. Design all primers across full alignment ---
  prof <- consensusProfile(aln_msa)
  oligos <- do.call(designOligos, c(list(prof), design_args))
  assays <- designAssays(oligos, tmDifferencePrimers = 5, length = c(min_amplicon, max_amplicon))
  assays_df <- as.data.frame(assays)
  assays_df$Assay <- seq_len(nrow(assays_df))
  
  all_results <- list()
  for (r in seq_len(nrow(gappy_regions))) {
    win_start <- gappy_regions$start[r]
    win_end   <- gappy_regions$end[r]
    gap_p     <- gappy_regions$gap_prop[r]
    # --- 3. Keep primer pairs whose amplicon covers the gappy window ---
    filtered <- assays_df %>%
      mutate(
        amplicon_start = pmin(startFwd, endRev),
        amplicon_end   = pmax(endFwd, startRev)
      ) %>%
      filter(amplicon_start <= win_start, amplicon_end >= win_end)
    if (nrow(filtered) == 0) next
    
    # --- 4. For >70% overlapping amplicons, keep only top 3 by score ---
    keep_idx <- c()
    used_regions <- list()
    for (i in seq_len(nrow(filtered))) {
      s1 <- filtered$amplicon_start[i]
      e1 <- filtered$amplicon_end[i]
      len1 <- e1 - s1 + 1
      overlap_counts <- 0
      for (j in seq_along(used_regions)) {
        s2 <- used_regions[[j]][1]
        e2 <- used_regions[[j]][2]
        overlap <- max(0, min(e1, e2) - max(s1, s2) + 1)
        overlap_frac <- overlap / min(len1, e2 - s2 + 1)
        if (overlap_frac > overlap_fraction) overlap_counts <- overlap_counts + 1
      }
      if (overlap_counts < max_per_amp_region) {
        keep_idx <- c(keep_idx, i)
        used_regions[[length(used_regions) + 1]] <- c(s1, e1)
      }
    }
    region_out <- filtered[keep_idx, ]
    if (nrow(region_out) > 0) {
      region_out$target_window_start <- win_start
      region_out$target_window_end <- win_end
      region_out$target_gap_prop <- gap_p
      all_results[[length(all_results) + 1]] <- region_out
    }
  }
  final_df <- bind_rows(all_results)
  if (nrow(final_df) == 0) {
    message("No valid amplicons found for the gappy regions.")
    return(final_df)
  }
  
  # --- 5. Compute predicted product length for each sample ---
  aln_mat <- as.matrix(aln_msa)
  seq_names <- rownames(aln_mat)
  amplicon_lengths <- do.call(rbind, lapply(seq_len(nrow(final_df)), function(i) {
    idx1 <- final_df$amplicon_start[i]
    idx2 <- final_df$amplicon_end[i]
    lens <- apply(aln_mat[, idx1:idx2, drop=FALSE], 1, function(x) sum(x != "-"))
    data.frame(
      Assay = final_df$Assay[i],
      Sequence = seq_names,
      AmpliconLength = lens
    )
  }))
  
  amplicon_lengths <- amplicon_lengths %>%
    group_by(Assay, Sequence) %>%
    summarise(AmpliconLength = first(AmpliconLength), .groups = "drop")
  amplicon_wide <- amplicon_lengths %>%
    pivot_wider(id_cols = Assay, names_from = Sequence, values_from = AmpliconLength)
  
  # --- 6. Output: Only keep key columns and amplicon lengths ---
  summary_df <- final_df %>%
    select(
      Assay,
      iupacSequenceFwd, startFwd, endFwd, tmMeanFwd, gcContentMeanFwd,
      iupacSequenceRev, startRev, endRev, tmMeanRev, gcContentMeanRev,
      amplicon_start, amplicon_end,
      target_window_start, target_window_end, target_gap_prop, score
    )
  out_df <- left_join(summary_df, amplicon_wide, by = "Assay") %>%
    arrange(target_window_start, score) %>%
    slice_head(n = output_top)
  return(out_df)
}

# ==== Example usage ====
design_args = list(
  lengthPrimer = c(18, 24),
  maxGapFrequency = 0.3,
  maxDegeneracyPrimer = 8,
  gcPrimer = c(0.40, 0.65),
  tmPrimer = c(45, 65),
  probe = FALSE
)

result <- findAndDesignGappyAmplicons("fish.fasta", design_args = design_args)


# this is good!
findAndDesignDivergentAmplicons <- function(
    fasta_file,
    window_size = 100,
    window_step = 100,
    top_n_windows = 3,
    min_amplicon = 50, max_amplicon = 600,
    overlap_fraction = 0.7, max_per_amp_region = 3,
    design_args = list(),
    output_top = 20,
    save_distance_csv = FALSE, # TRUE可将矩阵输出到当前目录
    prefix = "divergent"
) {
  # 1. Read and align
  raw_dna <- readDNAStringSet(fasta_file)
  alignment <- msa(raw_dna, method = "ClustalW")
  aln_msa <- DNAMultipleAlignment(alignment)
  aln_mat <- as.matrix(aln_msa)
  aln_len <- ncol(aln_mat)
  seq_names <- rownames(aln_mat)
  
  # --- 1. Sliding window: measure divergence ---
  window_info <- data.frame()
  dist_matrices <- list()
  for (start in seq(1, aln_len - window_size + 1, by = window_step)) {
    end <- start + window_size - 1
    win <- aln_mat[, start:end, drop=FALSE]
    seqs <- apply(win, 1, paste0, collapse = "")
    # 计算pairwise hamming距离
    dists <- stringdistmatrix(seqs, seqs, method = "hamming")
    rownames(dists) <- seq_names
    colnames(dists) <- seq_names
    mean_dist <- mean(dists[upper.tri(dists)])   # 平均距离
    max_dist  <- max(dists[upper.tri(dists)])    # 最大距离
    window_info <- rbind(window_info,
                         data.frame(start = start, end = end, size = window_size, mean_dist = mean_dist, max_dist = max_dist)
    )
    dist_matrices[[paste0(start, "_", end)]] <- dists
  }
  divergent_regions <- window_info %>%
    arrange(desc(mean_dist)) %>%
    slice_head(n = top_n_windows)
  
  if (nrow(divergent_regions) == 0) {
    message("No divergent regions found in alignment.")
    return(NULL)
  }
  
  # 可选：保存每个窗口的距离矩阵到csv
  if(save_distance_csv) {
    for(i in seq_len(nrow(divergent_regions))) {
      win_start <- divergent_regions$start[i]
      win_end   <- divergent_regions$end[i]
      key <- paste0(win_start, "_", win_end)
      outfn <- sprintf("%s_window_%s_dist.csv", prefix, key)
      write.csv(as.matrix(dist_matrices[[key]]), file = outfn)
    }
  }
  
  # --- 2. Design all primers across full alignment ---
  prof <- consensusProfile(aln_msa)
  oligos <- do.call(designOligos, c(list(prof), design_args))
  assays <- designAssays(oligos, tmDifferencePrimers = 5, length = c(min_amplicon, max_amplicon))
  assays_df <- as.data.frame(assays)
  assays_df$Assay <- seq_len(nrow(assays_df))
  
  all_results <- list()
  for (r in seq_len(nrow(divergent_regions))) {
    win_start <- divergent_regions$start[r]
    win_end   <- divergent_regions$end[r]
    mean_dist <- divergent_regions$mean_dist[r]
    max_dist  <- divergent_regions$max_dist[r]
    # --- 3. 保留扩增区覆盖该窗口的所有引物对 ---
    filtered <- assays_df %>%
      mutate(
        amplicon_start = pmin(startFwd, endRev),
        amplicon_end   = pmax(endFwd, startRev)
      ) %>%
      filter(amplicon_start <= win_start, amplicon_end >= win_end)
    if (nrow(filtered) == 0) next
    
    # --- 4. 对于>70%重叠的扩增子，每组只保留3个（按score优先） ---
    keep_idx <- c()
    used_regions <- list()
    for (i in seq_len(nrow(filtered))) {
      s1 <- filtered$amplicon_start[i]
      e1 <- filtered$amplicon_end[i]
      len1 <- e1 - s1 + 1
      overlap_counts <- 0
      for (j in seq_along(used_regions)) {
        s2 <- used_regions[[j]][1]
        e2 <- used_regions[[j]][2]
        overlap <- max(0, min(e1, e2) - max(s1, s2) + 1)
        overlap_frac <- overlap / min(len1, e2 - s2 + 1)
        if (overlap_frac > overlap_fraction) overlap_counts <- overlap_counts + 1
      }
      if (overlap_counts < max_per_amp_region) {
        keep_idx <- c(keep_idx, i)
        used_regions[[length(used_regions) + 1]] <- c(s1, e1)
      }
    }
    region_out <- filtered[keep_idx, ]
    if (nrow(region_out) > 0) {
      region_out$target_window_start <- win_start
      region_out$target_window_end <- win_end
      region_out$window_mean_dist <- mean_dist
      region_out$window_max_dist <- max_dist
      region_out$dist_matrix_key <- paste0(win_start, "_", win_end)
      all_results[[length(all_results) + 1]] <- region_out
    }
  }
  final_df <- bind_rows(all_results)
  if (nrow(final_df) == 0) {
    message("No valid amplicons found for the divergent regions.")
    return(list(table = NULL, matrices = dist_matrices))
  }
  
  # --- 5. 每个物种的预测扩增产物长度 ---
  amplicon_seqs <- do.call(rbind, lapply(seq_len(nrow(final_df)), function(i) {
    idx1 <- final_df$amplicon_start[i]
    idx2 <- final_df$amplicon_end[i]
    seqs <- apply(aln_mat[, idx1:idx2, drop=FALSE], 1, function(x) paste0(x, collapse=""))
    data.frame(
      Assay = final_df$Assay[i],
      Sequence = seq_names,
      AmpliconSeq = seqs
    )
  }))
  amplicon_seqs <- amplicon_seqs %>%
    group_by(Assay, Sequence) %>%
    summarise(AmpliconSeq = first(AmpliconSeq), .groups = "drop")
  amplicon_wide <- amplicon_seqs %>%
    pivot_wider(id_cols = Assay, names_from = Sequence, values_from = AmpliconSeq)
  
  
  # --- 6. 为每个Assay计算产物区域的pairwise距离矩阵 ---
  assay_dist_matrices <- list()
  for (i in seq_len(nrow(final_df))) {
    idx1 <- final_df$amplicon_start[i]
    idx2 <- final_df$amplicon_end[i]
    amplicon_region <- aln_mat[, idx1:idx2, drop=FALSE]
    amplicon_seqs <- apply(amplicon_region, 1, paste0, collapse = "")
    dists <- stringdistmatrix(amplicon_seqs, amplicon_seqs, method = "hamming")
    rownames(dists) <- seq_names
    colnames(dists) <- seq_names
    # 可选：按Assay编号或产物位置命名
    key <- paste0("Assay", final_df$Assay[i], "_", idx1, "_", idx2)
    assay_dist_matrices[[key]] <- dists
  }
  
  # --- 6. 精简输出 ---
  summary_df <- final_df %>%
    select(
      Assay,
      iupacSequenceFwd, startFwd, endFwd, tmMeanFwd, gcContentMeanFwd,
      iupacSequenceRev, startRev, endRev, tmMeanRev, gcContentMeanRev,
      amplicon_start, amplicon_end,
      target_window_start, target_window_end, window_mean_dist, window_max_dist, dist_matrix_key, score
    )
  out_df <- left_join(summary_df, amplicon_wide, by = "Assay") %>%
    arrange(target_window_start, score) %>%
    slice_head(n = output_top)
  return(list(table = out_df, assay_dist_matrices = assay_dist_matrices))
}

# ==== Example usage ====
design_args = list(
  lengthPrimer = c(18, 24),
  maxGapFrequency = 0.3,
  maxDegeneracyPrimer = 8,
  gcPrimer = c(0.40, 0.7),
  tmPrimer = c(45, 65),
  probe = FALSE
)

result <- findAndDesignDivergentAmplicons("./0618_fna/6432.FNA", design_args = design_args, save_distance_csv = TRUE)
if(!is.null(result$table)) {
  # 查看第1对引物的产物距离矩阵
  mat_name <- names(result$assay_dist_matrices)[1]
  print(mat_name)
  print(result$assay_dist_matrices[[mat_name]])
}

buildTreeForAmpliconRow <- function(result_df, row_id = 1, output_prefix = "assay_tree", bootstrap = 100) {
  meta_cols <- c(
    "Assay", "iupacSequenceFwd", "startFwd", "endFwd", "tmMeanFwd", "gcContentMeanFwd",
    "iupacSequenceRev", "startRev", "endRev", "tmMeanRev", "gcContentMeanRev",
    "amplicon_start", "amplicon_end", "target_window_start", "target_window_end",
    "window_mean_dist", "window_max_dist", "score", "dist_matrix_key"
  )
  sample_cols <- setdiff(colnames(result_df), meta_cols)
  seqs <- unlist(result_df[row_id, sample_cols])
  names(seqs) <- sample_cols

  seqs <- seqs[!is.na(seqs) & seqs != "" & grepl("^[ACGTNacgtn-]+$", seqs)]
  if (length(seqs) < 3) stop("有效序列少于3个，无法建树！")

  amplicon_set <- DNAStringSet(seqs)
  names(amplicon_set) <- names(seqs)
  fasta_file <- sprintf("%s_row%s.fasta", output_prefix, row_id)
  writeXStringSet(amplicon_set, fasta_file)

  # 用ape::read.dna + dist.dna建NJ树
  seqs2 <- read.dna(fasta_file, format = "fasta")
  dist_matrix <- dist.dna(seqs2, model = "K80")
  tree_nj <- nj(dist_matrix)

  # 可选：直接画NJ树
  pdf_file1 <- sprintf("%s_row%s_NJ.pdf", output_prefix, row_id)
  pdf(pdf_file1, width = 8, height = 8)
  plot(tree_nj, main = "Phylogenetic tree (NJ)")
  dev.off()

  # 若还要ML+bootstrap，则用phangorn::phyDat + ML优化
  phydat <- phyDat(seqs2, type = "DNA")
  fit <- pml(tree_nj, phydat)
  fit_opt <- optim.pml(fit, model = "HKY", rearrangement = "stochastic")
  bs <- bootstrap.pml(fit_opt, bs = bootstrap, optNni = TRUE)
  tree_file <- sprintf("%s_row%s_ML.nwk", output_prefix, row_id)
  ape::write.tree(fit_opt$tree, tree_file)
  pdf_file2 <- sprintf("%s_row%s_ML.pdf", output_prefix, row_id)
  pdf(pdf_file2, width = 8, height = 8)
  plotBS(midpoint(fit_opt$tree), bs, p = 50, main = sprintf("Assay row %s", row_id))
  dev.off()

  message("输出文件：\n", fasta_file, "\n", tree_file, "\n", pdf_file1, "\n", pdf_file2)
  invisible(list(fasta = fasta_file, nwk = tree_file, NJ_pdf = pdf_file1, ML_pdf = pdf_file2, tree = fit_opt$tree, bs = bs))
}

mytable <- result$table
buildTreeForAmpliconRow(mytable, row_id = 1, output_prefix = "amplicon", bootstrap = 100)

plotGelMapForAssayRow <- function(result_df, row_id = 1,
                                  marker = c(100, 200, 300, 400, 500, 600, 800, 1200),
                                  lane_names = NULL, band_col = "black", marker_col = "blue",
                                  band_width = 0.7) {
  meta_cols <- c(
    "Assay","iupacSequenceFwd","startFwd","endFwd","tmMeanFwd","gcContentMeanFwd",
    "iupacSequenceRev","startRev","endRev","tmMeanRev","gcContentMeanRev",
    "amplicon_start","amplicon_end","target_window_start","target_window_end",
    "target_gap_prop","score"
  )
  band_df <- result_df[row_id, setdiff(colnames(result_df), meta_cols)]
  band_tbl <- tibble(Lane = names(band_df), Size = as.numeric(band_df[1,]))
  band_tbl <- band_tbl[!is.na(band_tbl$Size) & band_tbl$Size > 0, ]
  if(is.null(lane_names)) lane_names <- band_tbl$Lane
  marker_tbl <- tibble(Lane = "Marker", Size = marker)
  band_tbl$Lane <- factor(band_tbl$Lane, levels = lane_names)
  marker_tbl$Lane <- factor(marker_tbl$Lane, levels = c("Marker"))
  all_lanes <- c("Marker", as.character(lane_names))
  band_tbl$Lane <- factor(band_tbl$Lane, levels = all_lanes)

  # 设定胶高范围（对数映射），大分子在上，小分子在下
  min_bp <- min(c(marker, band_tbl$Size))
  max_bp <- max(c(marker, band_tbl$Size))
  gel_top <- log10(max_bp)   # 胶上缘（最小迁移，最大bp）
  gel_bottom <- log10(min_bp) # 胶下缘（最大迁移，最小bp）

  # band height (像素/相对值)
  bar_height <- 0.04
  # 映射y轴为-log10(bp)模拟迁移距离
  band_tbl$y <- -log10(band_tbl$Size)
  marker_tbl$y <- -log10(marker_tbl$Size)

  # ggplot主图
  p <- ggplot() +
    geom_tile(data = marker_tbl, aes(x = Lane, y = y, width = band_width, height = bar_height),
              fill = marker_col, alpha = 0.7) +
    geom_tile(data = band_tbl, aes(x = Lane, y = y, width = band_width, height = bar_height),
              fill = band_col, alpha = 0.9) +
    # 刻度线和y轴反向
    scale_y_continuous(
      breaks = -log10(marker),
      labels = marker,
      trans = "identity",
      limits = c(-log10(max_bp)-0.05, -log10(min_bp)+0.05)
    ) +
    labs(y = "Fragment size (bp)", x = "Lane", row_id) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank())
  print(p)
  invisible(p)
}

plotGelMapForAssayRow(result, row_id = 1)

# # -------------
# 
# g_name <- ""
# 
# extract_shortest_interval_between_genes(g_name, input_dir = "./", output_fasta = paste0(g_name, ".fasta"))
# 
# align_and_calc_distance(paste0(g_name,".fasta"), paste0(g_name, "_distance.tsv"))
# 
# result <- classify_species_by_distance(paste0(g_name, "_distance.tsv"))
# 
# write.table(result, paste0(g_name, "_classification.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# 
# # -------------
# 
# g_name1 <- ""
# g_name2 <- ""
# 
# extract_shortest_interval_between_genes(g_name1, g_name2, input_dir = "./", output_fasta = paste0(g_name1, "_", g_name2, ".fasta"))
# 
# align_and_calc_distance(paste0(g_name1, "_", g_name2, ".fasta"), paste0(g_name1, "_", g_name2,  "_distance.tsv"))
# 
# result <- classify_species_by_distance(paste0(g_name1, "_", g_name2,  "_distance.tsv"))
# 
# write.table(result, paste0(g_name1, "_", g_name2,  "_classification.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
# 
# # -------------
# 
# results <- combine_and_evaluate_classification("./")
# write.table(results, file = "gene_combination_summary.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# 
# 

# # -------------

# process_all_common_regions(input_dir = "./0526test_split", mode = "length_indel")

# --------

process_all_common_FNA <- function(input_dir = ".", output_dir = ".", t_offhit = 0.6, r_offhit = 0.7, t_indist = 0.01, mode = c("distance", "length_indel")) {
  mode <- match.arg(mode)
  
  # 1️⃣ 直接获取目录下所有 .fna 文件
  fna_files <- list.files(input_dir, pattern = "\\.FNA$", full.names = TRUE)
  if (length(fna_files) == 0) {
    stop("未找到任何 .fna 文件。请检查 input_dir 路径。")
  }
  
  # 2️⃣ 用 .fna 文件名（不含后缀）作为候选区域列表
  all_targets <- basename(fna_files)
  all_targets <- sub("\\.FNA$", "", all_targets)
  
  for (i in seq_along(fna_files)) {
    tag <- all_targets[i]
    fasta_file <- fna_files[i]
    dist_file <- file.path(output_dir, paste0(tag, "_distance.tsv"))
    class_file <- file.path(output_dir, paste0(tag, "_classification.tsv"))
    
    # 不再提取基因或区间，直接用 fna 文件做后续分析
    if (mode == "distance") {
      align_and_calc_distance(fasta_file, dist_file)
      result <- classify_species_by_distance(dist_file,
                                             offhit_threshold = t_offhit,
                                             offhit_ratio = r_offhit,
                                             indist_threshold = t_indist)
      write.table(result, file = class_file, sep = "\t", row.names = FALSE, quote = FALSE)
    } else if (mode == "length_indel") {
      indel_result <- find_visible_indel_windows(fasta_file)
      if (!is.null(indel_result)) {
        write.table(indel_result, file = file.path(output_dir, paste0(tag, "_visibleindel.tsv")),
                    sep = "\t", row.names = FALSE, quote = FALSE)
        screen_fasta_for_pcr_indels(fasta_file, output_dir = output_dir)
      }
    }
  }
  
  # 3️⃣ 汇总/结束提示
  if (mode == "distance") {
    summary <- combine_and_evaluate_classification(output_dir)
    write.table(summary, file = file.path(output_dir, "combined_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    cat("🎉 所有区域已处理并评价。\n")
  } else {
    cat("🎉 所有区域已扫描 indel。\n")
  }
}

process_all_common_FNA(input_dir = "./0618_fna", mode = "distance")

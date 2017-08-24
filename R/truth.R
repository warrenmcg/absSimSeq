# actual function to compare simulation results to the ground truth
sleuth_sim_to_truth <- function(index, sleuth_dir = ".", gene_mode = FALSE,
                         group_ids = NULL, final_results = NULL,
                         prefixes = NULL, test = "wt") {
  stopifnot(test %in% c("wt", "lrt"))
  message("comparing run #", index)
  if (is.null(prefixes)) {
    out_dir <- sleuth_dir
    in_dir <- file.path(sleuth_dir, paste0("sleuth_out_run", index))
    prefix <- file.path(in_dir, paste0("run", index))
    if (gene_mode) {
      trans_file <- paste0(prefix, "Gene.RData")
      alr_file <- paste0(prefix, "_alrGene.RData")
    } else {
      trans_file <- paste0(prefix, ".RData")
      alr_file <- paste0(prefix, "_alr.RData")
    }
  } else {
    trans_file <- file.path(sleuth_dir, paste0(prefixes[1], ".RData"))
    alr_file <- file.path(sleuth_dir, paste0(prefixes[2], ".RData"))
    prefix <- prefixes[3]
    out_dir <- dirname(prefixes[3])
  }
  
  message("loading ", trans_file)
  load(trans_file)
  if (test == "wt") {
    old_results <- sleuth::sleuth_results(sleuth.obj, "conditionexperiment")
  } else {
    old_results <- sleuth::sleuth_results(sleuth.obj, "reduced:full",
                                          test_type = "lrt")
  }
  if (gene_mode) {
    old_counts <- sleuth:::spread_abundance_by(sleuth.obj$obs_norm,
                                               "scaled_reads_per_base")
  } else {
    old_counts <- sleuth:::spread_abundance_by(sleuth.obj$obs_norm,
                                               "est_counts")
  }
  rm(sleuth.obj)
  old_ctr_mean <- apply(old_counts[, which(group_ids==1)],
                        1, mean)
  old_exp_mean <- apply(old_counts[, which(group_ids==2)],
                        1, mean)
  old_means <- data.frame(ctr_mean = old_ctr_mean, exp_mean = old_exp_mean)
  rownames(old_means) <- rownames(old_counts)
  old_results <- merge(old_results, old_means,
                       by.x = "target_id", by.y = "row.names", all = T)
  message("loading ", alr_file)
  load(alr_file)
  if (test == "wt") {
    alr_results <- sleuth::sleuth_results(sleuth.obj, "conditionexperiment")
  } else {
    alr_results <- sleuth::sleuth_results(sleuth.obj, "reduced:full",
                                          test_type = "lrt")
  }
  alr_tpms <- sleuth.obj$bs_summary$obs_tpm
  rm(sleuth.obj)
  alr_ctr_mean <- apply(alr_tpms[, which(group_ids==1)],
                        1, mean)
  alr_exp_mean <- apply(alr_tpms[, which(group_ids==2)],
                        1, mean)
  alr_means <- data.frame(ctr_mean = alr_ctr_mean, exp_mean = alr_exp_mean)
  rownames(alr_means) <- rownames(alr_tpms)
  alr_results <- merge(alr_results, alr_means,
                       by.x = "target_id", by.y = "row.names", all = T)
  alr_result <- final_results$alr_data[[index]]
  alr_result <- as.data.frame(alr_result)
  names(alr_result)[1:2] <- c("alr_ctr_ratio", "alr_exp_ratio")
  #rownames(alr_result) <- substr(rownames(alr_result), 1, 15)
  message("comparing these results with the truth")
  main_result <- final_results$results[[index]]
  if (gene_mode) {
    main_result <- main_result[!is.na(main_result[, 1]), ]
    row.names(main_result) <- main_result[, 1]
    main_result$ensembl_gene_id <- NULL
    names(main_result)[1:2] <- c("ctr_copy_numbers", "exp_copy_numbers")
  } else {
    main_result <- as.data.frame(main_result)
    names(main_result)[4:9] <- c("ctr_copy_numbers", "exp_copy_numbers",
                                 "ctr_tpms", "exp_tpms",
                                 "ctr_reads", "exp_reads")
  }
  if (test == "wt") {
    cols2include <- c("target_id", "gene_symbol", "pval", "qval", "b", "se_b",
                      "ctr_mean", "exp_mean")
  } else {
    cols2include <- c("target_id", "gene_symbol", "pval", "qval", "test_stat", "rss",
                      "ctr_mean", "exp_mean")
  }
  old_comparison <- merge(old_results[, cols2include], main_result,
                          by.x = "target_id", by.y = "row.names",
                          all = T)
  if(any(is.na(old_comparison$qval) &
         old_comparison$ctr_copy_numbers == 0)) {
    old_comparison <- old_comparison[-which(
      is.na(old_comparison$qval) &
        old_comparison$ctr_copy_numbers == 0), ]
  }
  if (test == "wt") {
    old_comparison <- old_comparison[order(old_comparison$qval,
                                           old_comparison$pval,
                                           -old_comparison$b), ]
  } else {
    old_comparison <- old_comparison[order(old_comparison$qval,
                                           old_comparison$pval,
                                           old_comparison$rss), ]
  }
  old_comparison$quartile <- cut(old_comparison$ctr_copy_numbers, 
                                 summary(old_comparison$ctr_copy_numbers)[-4])
  if (test == "wt") {
    old_comparison$de_status <- ifelse(old_comparison$abs_fold_change == 1,
                                       0, ifelse(
                                         old_comparison$abs_fold_change >= 1,
                                         1, -1))
    old_comparison$direction <- ifelse(is.na(old_comparison$qval), 0,
                            ifelse(old_comparison$b > 0, 1, -1))
  } else {
    old_comparison$de_status <- ifelse(old_comparison$abs_fold_change == 1,
                                       0, 1)
    old_comparison$direction <- ifelse(is.na(old_comparison$qval), 0, 1)
  }
  old_comparison$decision <- ifelse(old_comparison$ctr_copy_numbers == 0 |
                                      old_comparison$de_status == 0,
                                    0, ifelse(old_comparison$direction != 0 &
                                                old_comparison$de_status ==
                                                old_comparison$direction,
                                              1, 0))
  old_comparison$tpr <- old_comparison$fpr <- old_comparison$fdr <- NA
  old_comparison$tpr <- cumsum(old_comparison$decision) / 
    sum(old_comparison$decision, na.rm=T)
  old_comparison$fpr <- cumsum(
    -1*(old_comparison$decision-1)) /
    (-1*(sum(old_comparison$decision-1, na.rm=T)))
  old_comparison$fdr <- cumsum(-1*(old_comparison$decision-1)) / 
    (seq(length(old_comparison$decision)))
  
  alr_comparison <- merge(alr_results[, cols2include], alr_result,
                          by.x = "target_id", by.y = "row.names",
                          all = T)
  alr_comparison <- merge(alr_comparison, main_result,
                            by.x = "target_id", by.y = "row.names")
  # exclude any targets with zero true copy numbers and was filtered out
  if (any(is.na(alr_comparison$qval) &
          alr_comparison$ctr_copy_numbers == 0)) {
    alr_comparison <- alr_comparison[-which(
      is.na(alr_comparison$qval) &
        alr_comparison$ctr_copy_numbers == 0), ]
  }
  # targets with zero true copy numbers but had estimated expression are
  # considered false positives (alr_diff == 0)
  alr_comparison$alr_diff[is.na(alr_comparison$alr_diff)] <- 0
  if (test == "wt") {
    alr_comparison <- alr_comparison[order(alr_comparison$qval,
                                           alr_comparison$pval,
                                           -alr_comparison$b), ]
  } else {
    alr_comparison <- alr_comparison[order(alr_comparison$qval,
                                           alr_comparison$pval,
                                           alr_comparison$rss), ]
  }
  alr_comparison$quartile <- cut(alr_comparison$ctr_copy_numbers, 
                                 summary(alr_comparison$ctr_copy_numbers)[-4])
  if (test == "wt") {
    alr_comparison$de_status <- ifelse(abs(alr_comparison$alr_diff) < 1e-5,
                                       0, ifelse(
                                         alr_comparison$alr_diff > 0,
                                         1, -1))
    alr_comparison$direction <- ifelse(is.na(alr_comparison$qval), 0,
                            alr_comparison$b > 0, 1, -1))
  } else {
    alr_comparison$de_status <- ifelse(abs(alr_comparison$alr_diff) < 1e-5,
                                       0, 1)
    alr_comparison$direction <- ifelse(is.na(alr_comparison$qval), 0, 1)
  }
  alr_comparison$decision <- ifelse(alr_comparison$ctr_copy_numbers == 0 |
                                      alr_comparison$de_status == 0,
                                    0, ifelse(alr_comparison$direction != 0 &
                                                alr_comparison$de_status ==
                                                alr_comparison$direction,
                                              1, 0))
  alr_comparison$tpr <- alr_comparison$fpr <- alr_comparison$fdr <- NA
  alr_comparison$tpr <- cumsum(alr_comparison$decision) / 
    sum(alr_comparison$decision, na.rm=T)
  alr_comparison$fpr <- cumsum(-1*(alr_comparison$decision-1)) /
    (-1*(sum(alr_comparison$decision-1, na.rm=T)))
  alr_comparison$fdr <- cumsum(-1*(alr_comparison$decision-1)) / 
    (seq(length(alr_comparison$decision)))
  
  plot_data <- data.frame(method = "alr", TPR = alr_comparison$tpr,
                          FPR = alr_comparison$fpr, FDR = alr_comparison$fdr)
  plot_data2 <- data.frame(method = "sleuth", TPR = old_comparison$tpr,
                           FPR = old_comparison$fpr, FDR = old_comparison$fdr)
  plot_data <- rbind(plot_data, plot_data2)
  plot_data <- plot_data[which(!is.na(plot_data$TPR)), ]
  g <- ggplot2::ggplot(plot_data, ggplot2::aes(x = FPR, y = TPR,
                                               color = method)) + 
    ggplot2::geom_line() + 
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype="dashed",
                         alpha = 0.3) +
    ggplot2::theme(axis.ticks.length = ggplot2::unit(0.5,"cm"),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   strip.text.x = ggplot2::element_text(size = 16,
                                                        face = "bold"),
                   legend.text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_text(size = 20))
  if (gene_mode) {
    suffix <- "_geneRoc.pdf"
  } else {
    suffix <- "_roc.pdf"
  }
  
  prefix <- basename(prefix)
  pdf(file.path(out_dir, paste0(prefix, suffix)))
  print(g)
  invisible(dev.off())
  
  g2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = FDR, y = TPR,
                                                color = method)) + 
    ggplot2::geom_path(size = 0.8, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0.05, linetype="dashed",
                        alpha = 0.3) + 
    ggplot2::geom_vline(xintercept = 0.1, linetype="dashed",
                        alpha = 0.3) +
    ggplot2::theme(axis.ticks.length = ggplot2::unit(0.5,"cm"),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   strip.text.x = ggplot2::element_text(size = 16,
                                                        face = "bold"),
                   legend.text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_text(size = 20))
  if (gene_mode) {
    suffix <- "_geneSensVsFDR.pdf"
  } else {
    suffix <- "_sensVsFDR.pdf"
  }
  pdf(file.path(out_dir, paste0(prefix, suffix)))
  print(g2)
  invisible(dev.off())
  
  list(sleuth = old_comparison, alr = alr_comparison)
}

other_sim_to_truth <- function(index, in_dir = ".", gene_mode = FALSE,
                               group_ids = NULL, final_results = NULL,
                               prefixes = NULL, tool = "DESeq2") {
  message("comparing run #", index)
  s_prefix <- file.path(in_dir, prefixes[1])
  k_prefix <- file.path(in_dir, prefixes[2])
  prefix <- prefixes[3]
  out_dir <- dirname(prefix)

  if (tool == "ALDEx2we") {
    type <- "alr"
    s_file <- paste(s_prefix, "aldex2AllResults.txt", sep = "_")
    k_file <- paste(k_prefix, "aldex2AllResults.txt", sep = "_")
    cols2include <- c("target_id", "gene_symbol", "rab.win.control",
                      "rab.win.experiment", "effect", "we.ep", "we.eBH")
  } else if (tool == "ALDEx2overlap") {
    type <- "alr"
    s_file <- paste(s_prefix, "aldex2AllResults.txt", sep = "_")
    k_file <- paste(k_prefix, "aldex2AllResults.txt", sep = "_")
    cols2include <- c("target_id", "gene_symbol", "rab.win.control",
                      "rab.win.experiment", "we.ep", "we.eBH", "effect", "overlap")
  } else if (tool == "ALDEx2wi") {
    type <- "alr"
    s_file <- paste(s_prefix, "aldex2AllResults.txt", sep = "_")
    k_file <- paste(k_prefix, "aldex2AllResults.txt", sep = "_")
    cols2include <- c("target_id", "gene_symbol", "rab.win.control",
                      "rab.win.experiment", "effect", "wi.ep", "wi.eBH")
  } else if (tool == "DESeq2") {
    type <- "main"
    s_file <- paste(s_prefix, "DESeqAllResults.txt", sep = "_")
    k_file <- paste(k_prefix, "DESeqAllResults.txt", sep = "_")
    cols2include <- c("target_id", "gene_symbol", "baseMean",
                      "log2FoldChange", "pvalue", "padj")
  } else if (tool == "limma") {
    type <- "main"
    s_file <- paste(s_prefix, "limmaAllResults.txt", sep = "_")
    k_file <- paste(k_prefix, "limmaAllResults.txt", sep = "_")
    cols2include <- c("target_id", "gene_symbol", "AveExpr",
                      "logFC", "P.Value", "adj.P.Val")
  }

  message('loading salmon results')
  s_results <- read.table(s_file, sep = "\t", header = T, quote = "")
  message('loading kallisto results')
  k_results <- read.table(k_file, sep = "\t", header = T, quote = "")

  message('retrieving the truth')
  main_result <- final_results$results[[index]]
  main_result <- as.data.frame(main_result)
  if (gene_mode) {
    main_result <- main_result[!is.na(main_result[, 1]), ]
    row.names(main_result) <- main_result[, 1]
    main_result$ensembl_gene_id <- NULL
    names(main_result)[1:2] <- c("ctr_copy_numbers", "exp_copy_numbers")
  } else {
    names(main_result)[4:9] <- c("ctr_copy_numbers", "exp_copy_numbers",
                                 "ctr_tpms", "exp_tpms",
                                 "ctr_reads", "exp_reads")
  }
  main_result$target_id <- rownames(main_result)

  if (type == "alr") {
    alr_result <- final_results$alr_data[[index]]
    alr_result <- as.data.frame(alr_result)
    names(alr_result)[1:2] <- c("alr_ctr_ratio", "alr_exp_ratio")
    alr_result$target_id <- rownames(alr_result)
    main_result <- merge(data.table::as.data.table(alr_result),
                         data.table::as.data.table(main_result),
                         by = "target_id", all = T)
    main_result <- as.data.frame(main_result)
  }

  message('merging results with the truth')
  s_comparison <- merge(data.table::as.data.table(s_results[, cols2include]),
                        data.table::as.data.table(main_result),
                        by = "target_id", all = T)
  s_comparison <- as.data.frame(s_comparison)
  s_comparison$quartile <- cut(s_comparison$ctr_copy_numbers,
                               unique(summary(s_comparison$ctr_copy_numbers)[-4]))
  k_comparison <- merge(data.table::as.data.table(k_results[, cols2include]),
                        data.table::as.data.table(main_result),
                        by = "target_id", all = T)
  k_comparison <- as.data.frame(k_comparison)
  k_comparison$quartile <- cut(k_comparison$ctr_copy_numbers,
                               unique(summary(
                                 k_comparison$ctr_copy_numbers)[-4]))

  ## order the comparisons by p-value
  if (tool == "ALDEx2we") {
    s_comparison <- s_comparison[order(s_comparison$we.eBH,
                                       s_comparison$we.ep,
                                       -s_comparison$effect), ]
    k_comparison <- k_comparison[order(k_comparison$we.eBH,
                                       k_comparison$we.ep,
                                       -k_comparison$effect), ]
    if(any(is.na(s_comparison$we.eBH) &
           s_comparison$ctr_copy_numbers == 0)) {
      s_comparison <- s_comparison[-which(
        is.na(s_comparison$we.eBH) &
          s_comparison$ctr_copy_numbers == 0), ]
    }
    if(any(is.na(k_comparison$we.eBH) &
           k_comparison$ctr_copy_numbers == 0)) {
      k_comparison <- k_comparison[-which(
        is.na(k_comparison$we.eBH) &
          k_comparison$ctr_copy_numbers == 0), ]
    }
  } else if (tool == "ALDEx2overlap") {
    s_comparison <- s_comparison[order(s_comparison$overlap,
                                       s_comparison$we.ep,
                                       -s_comparison$effect), ]
    k_comparison <- k_comparison[order(k_comparison$overlap,
                                       k_comparison$we.ep,
                                       -k_comparison$effect), ]
    if(any(is.na(s_comparison$overlap) &
           s_comparison$ctr_copy_numbers == 0)) {
      s_comparison <- s_comparison[-which(
        is.na(s_comparison$overlap) &
          s_comparison$ctr_copy_numbers == 0), ]
    }
    if(any(is.na(k_comparison$overlap) &
           k_comparison$ctr_copy_numbers == 0)) {
      k_comparison <- k_comparison[-which(
        is.na(k_comparison$overlap) &
          k_comparison$ctr_copy_numbers == 0), ]
    }
  } else if (tool == "ALDEx2wi") {
    s_comparison <- s_comparison[order(s_comparison$wi.eBH,
                                       s_comparison$wi.ep,
                                       -s_comparison$effect), ]
    k_comparison <- k_comparison[order(k_comparison$wi.eBH,
                                       k_comparison$wi.ep,
                                       -k_comparison$effect), ]
    if(any(is.na(s_comparison$wi.eBH) &
           s_comparison$ctr_copy_numbers == 0)) {
      s_comparison <- s_comparison[-which(
        is.na(s_comparison$wi.eBH) &
          s_comparison$ctr_copy_numbers == 0), ]
    }
    if(any(is.na(k_comparison$wi.eBH) &
           k_comparison$ctr_copy_numbers == 0)) {
      k_comparison <- k_comparison[-which(
        is.na(k_comparison$wi.eBH) &
          k_comparison$ctr_copy_numbers == 0), ]
    }
  } else if (tool == "DESeq2") {
    s_comparison <- s_comparison[order(s_comparison$padj,
                                       s_comparison$pvalue,
                                       -s_comparison$log2FoldChange), ]
    k_comparison <- k_comparison[order(k_comparison$padj,
                                       k_comparison$pvalue,
                                       -k_comparison$log2FoldChange), ]
    if(any(is.na(s_comparison$padj) &
           s_comparison$ctr_copy_numbers == 0)) {
      s_comparison <- s_comparison[-which(
        is.na(s_comparison$padj) &
          s_comparison$ctr_copy_numbers == 0), ]
    }
    if(any(is.na(k_comparison$padj) &
           k_comparison$ctr_copy_numbers == 0)) {
      k_comparison <- k_comparison[-which(
        is.na(k_comparison$padj) &
          k_comparison$ctr_copy_numbers == 0), ]
    }
    dir_col <- "log2FoldChange"
  } else if (tool == "limma") {
    s_comparison <- s_comparison[order(s_comparison$adj.P.Val,
                                       s_comparison$P.Value,
                                       -s_comparison$logFC), ]
    k_comparison <- k_comparison[order(k_comparison$adj.P.Val,
                                       k_comparison$P.Value,
                                       -k_comparison$logFC), ]
    if(any(is.na(s_comparison$adj.P.Val) &
           s_comparison$ctr_copy_numbers == 0)) {
      s_comparison <- s_comparison[-which(
        is.na(s_comparison$adj.P.Val) &
          s_comparison$ctr_copy_numbers == 0), ]
    }
    if(any(is.na(k_comparison$adj.P.Val) &
           k_comparison$ctr_copy_numbers == 0)) {
      k_comparison <- k_comparison[-which(
        is.na(k_comparison$adj.P.Val) &
          k_comparison$ctr_copy_numbers == 0), ]
    }
    dir_col <- "logFC"
  }

  if (type == "alr") {
    s_comparison$alr_diff[is.na(s_comparison$alr_diff)] <- 0
    k_comparison$alr_diff[is.na(k_comparison$alr_diff)] <- 0
    num_col <- "alr_diff"
    dir_col <- "effect"
  } else {
    num_col <- "abs_fold_change"
  }

  ## Common code for all tools to process the decision of True vs False Positive
  message('comparing salmon results to truth')
  s_comparison$de_status <- ifelse(s_comparison[, num_col] == 1,
                                   0, ifelse(
                                     s_comparison[, num_col] >= 1,
                                     1, -1))
  s_comparison$direction <- ifelse(is.na(s_comparison[, dir_col]), 0,
                                   ifelse(s_comparison[, dir_col] > 0, 1, -1))
  s_comparison$decision <- ifelse(s_comparison$ctr_copy_numbers == 0 |
                                    s_comparison$de_status == 0,
                                  0, ifelse(s_comparison$direction != 0 &
                                              s_comparison$de_status ==
                                              s_comparison$direction,
                                            1, 0))
  s_comparison$tpr <- s_comparison$fpr <- s_comparison$fdr <- NA
  s_comparison$tpr <- cumsum(s_comparison$decision) /
    sum(s_comparison$decision, na.rm=T)
  s_comparison$fpr <- cumsum(-1*(s_comparison$decision-1)) /
    (-1*(sum(s_comparison$decision-1, na.rm=T)))
  s_comparison$fdr <- cumsum(-1*(s_comparison$decision-1)) /
    (seq(length(s_comparison$decision)))

  message('comparing kallisto results to truth')
  k_comparison$de_status <- ifelse(k_comparison[, num_col] == 1,
                                   0, ifelse(
                                     k_comparison[, num_col] >= 1,
                                     1, -1))
  k_comparison$direction <- ifelse(is.na(k_comparison[, dir_col]), 0,
                                   ifelse(k_comparison[, dir_col] > 0, 1, -1))
  k_comparison$decision <- ifelse(k_comparison$ctr_copy_numbers == 0 |
                                    k_comparison$de_status == 0,
                                  0, ifelse(k_comparison$direction != 0 &
                                              k_comparison$de_status ==
                                              k_comparison$direction,
                                            1, 0))
  k_comparison$tpr <- k_comparison$fpr <- k_comparison$fdr <- NA
  k_comparison$tpr <- cumsum(k_comparison$decision) /
    sum(k_comparison$decision, na.rm=T)
  k_comparison$fpr <- cumsum(-1*(k_comparison$decision-1)) /
    (-1*(sum(k_comparison$decision-1, na.rm=T)))
  k_comparison$fdr <- cumsum(-1*(k_comparison$decision-1)) /
    (seq(length(k_comparison$decision)))

  message('plotting the results')
  plot_data <- data.frame(method = "salmon", TPR = s_comparison$tpr,
                          FPR = s_comparison$fpr, FDR = s_comparison$fdr)
  plot_data2 <- data.frame(method = "kallisto", TPR = k_comparison$tpr,
                           FPR = k_comparison$fpr, FDR = k_comparison$fdr)
  plot_data <- rbind(plot_data, plot_data2)
  plot_data <- plot_data[which(!is.na(plot_data$TPR)), ]
  g <- ggplot2::ggplot(plot_data, ggplot2::aes(x = FPR, y = TPR,
                                               color = method)) +
    ggplot2::geom_line() +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype="dashed",
                         alpha = 0.3) +
    ggplot2::theme(axis.ticks.length = ggplot2::unit(0.5,"cm"),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   strip.text.x = ggplot2::element_text(size = 16,
                                                        face = "bold"),
                   legend.text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_text(size = 20))
  if (gene_mode) {
    suffix <- "_geneRoc.pdf"
  } else {
    suffix <- "_roc.pdf"
  }

  prefix <- basename(prefix)
  prefix <- paste(prefix, tool, sep = "_")
  pdf(file.path(out_dir, paste0(prefix, suffix)))
  print(g)
  invisible(dev.off())

  g2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = FDR, y = TPR,
                                                color = method)) +
    ggplot2::geom_path(size = 0.8, alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0.05, linetype="dashed",
                        alpha = 0.3) +
    ggplot2::geom_vline(xintercept = 0.1, linetype="dashed",
                        alpha = 0.3) +
    ggplot2::geom_vline(xintercept = 0.1, linetype="dashed",
                        alpha = 0.3) +
    ggplot2::theme(axis.ticks.length = ggplot2::unit(0.5,"cm"),
                   axis.text.x = ggplot2::element_text(angle = 90,
                                                       vjust = 0.5,
                                                       size = 16),
                   axis.text.y = ggplot2::element_text(size = 16),
                   axis.title.x = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   axis.title.y = ggplot2::element_text(size = 20,
                                                        face = "bold"),
                   strip.text.x = ggplot2::element_text(size = 16,
                                                        face = "bold"),
                   legend.text = ggplot2::element_text(size = 16),
                   legend.title = ggplot2::element_text(size = 20))
  if (gene_mode) {
    suffix <- "_geneSensVsFDR.pdf"
  } else {
    suffix <- "_sensVsFDR.pdf"
  }
  pdf(file.path(out_dir, paste0(prefix, suffix)))
  print(g2)
  invisible(dev.off())

  list(salmon = s_comparison, kallisto = k_comparison)
}

compare_sim_to_truth <- function(final_results, in_dir = ".",
                                 de_probs, dir_probs,
                                 num_reps = c(10, 10), prefixes = NULL,
                                 gene_mode = FALSE, test = "wt",
                                 tool = "sleuth") {
  group_ids <- c(rep(1, num_reps[1]), rep(2, num_reps[2]))
  if (tool == "sleuth") {
    if (is.null(prefixes)) {
      sim2truth <- lapply(seq_along(final_results$results), sleuth_sim_to_truth,
                          sleuth_dir = in_dir, gene_mode = gene_mode,
                          group_ids = group_ids, final_results = final_results,
                          test = test)
    } else {
      sim2truth <- lapply(seq_along(final_results$results), function(x) {
        sleuth_sim_to_truth(index = x, sleuth_dir = in_dir,
                            gene_mode = gene_mode, group_ids = group_ids,
                            final_results = final_results,
                            prefixes = prefixes[[x]], test = test)
      })
    }
  } else {
    sim2truth <- lapply(seq_along(final_results$results), function(x) {
      other_sim_to_truth(index = x, sleuth_dir = in_dir, gene_mode = gene_mode,
                   group_ids = group_ids, final_results = final_results,
                   prefixes = prefixes[[x]], test = test, tool = tool)
    })
  }
  sim2truth
}

get_gene_truth <- function(final_results = NULL,
                           host = "dec2016.archive.ensembl.org",
                           species = "hsapiens",
                           denom = "ENSG00000111640") {
  test_ids <- head(rownames(final_results$results[[1]]$copy_numbers_per_cell))
  if (nchar(test_ids[1])>15) {
    message("loading annotations from biomaRt")
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                             host = host,
                             dataset = paste(species, "gene_ensembl", sep = "_"))
    annos <- biomaRt::getBM(mart = mart,
                            attributes = c("ensembl_transcript_id",
                                           "transcript_version",
                                           "ensembl_gene_id"))
    annos$target_id <- paste(annos$ensembl_transcript_id,
                             annos$transcript_version,
                             sep = ".")
    annos[, c("transcript_version", "ensembl_transcript_id")] <- NULL
  } else {
    message("loading annotations from biomaRt")
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                             host = host,
                             dataset = paste(species, "gene_ensembl", sep = "_"))
    annos <- biomaRt::getBM(mart = mart,
                            attributes = c("ensembl_transcript_id",
                                           "ensembl_gene_id"))
    annos$target_id <- annos$ensembl_transcript_id
    annos$ensembl_transcript_id <- NULL
  }
  gene_results <- lapply(seq(length(final_results$results)), function(index) {
    annos <- data.table::as.data.table(annos)
    result <- final_results$results[[index]]
    expected_reads <- result$expected_reads
    copy_nums <- as.data.frame(result$copy_numbers_per_cell)
    copy_nums$target_id <- row.names(copy_nums)
    copy_nums <- data.table::as.data.table(copy_nums)
    copy_nums <- merge(copy_nums, annos, by = "target_id", all.x = T)
    gene_copy_nums <- copy_nums[, .(ctr_copy_nums = sum(ctr_copy_numbers),
                                    exp_copy_nums = sum(exp_copy_numbers)),
                                by = "ensembl_gene_id"]
    gene_copy_nums[, abs_fold_changes := exp_copy_nums / ctr_copy_nums]
    gene_copy_nums[is.na(abs_fold_changes), abs_fold_changes := 1]
    gene_copy_nums[, ctr_tpms := ctr_copy_nums / sum(ctr_copy_nums) * 10^6]
    gene_copy_nums[, exp_tpms := exp_copy_nums / sum(exp_copy_nums) * 10^6]
    return(as.data.frame(gene_copy_nums))
  })
  
  alr_results <- lapply(gene_results, function(result) {
    result <- result[!is.na(result[, 1]), ]
    rownames(result) <- result[, 1]
    abs_copy_numbers <- result[, c(2, 3)]
    rel_tpms <- result[, c(5, 6)]
    calculate_rel_consistency(abs_copy_numbers, rel_tpms, denom)
  })
  return(list(results = gene_results, alr_data = alr_results))
}

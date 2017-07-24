# actual function to compare simulation results to the ground truth
sim_to_truth <- function(index, sleuth_dir = ".", gene_mode = FALSE,
                         group_ids = NULL, final_results = NULL,
                         prefixes = NULL) {
  message("comparing run #", index)
  if (is.null(prefixes)) {
    out_dir <- file.path(sleuth_dir, paste0("sleuth_out_run", index))
    prefix <- file.path(out_dir, paste0("run", index))
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
    prefix <- basename(prefixes[2])
  }
  
  message("loading ", trans_file)
  load(trans_file)
  old_results <- sleuth::sleuth_results(sleuth.obj, "conditionexperiment")
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
  alr_results <- sleuth::sleuth_results(sleuth.obj, "conditionexperiment")
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
  #rownames(main_result) <- substr(rownames(main_result), 1, 15)
  cols2include <- c("target_id", "gene_symbol", "pval", "qval", "b", "se_b",
                    "ctr_mean", "exp_mean")
  old_comparison <- merge(old_results[, cols2include], main_result,
                          by.x = "target_id", by.y = "row.names",
                          all = T)
  if(any(is.na(old_comparison$qval) &
         old_comparison$ctr_copy_numbers == 0)) {
    old_comparison <- old_comparison[-which(
      is.na(old_comparison$qval) &
        old_comparison$ctr_copy_numbers == 0), ]
  }
  old_comparison <- old_comparison[order(old_comparison$qval,
                                         old_comparison$pval,
                                         old_comparison$b), ]
  old_filtered <- which(!is.na(old_comparison$pval))
  old_comparison$quartile <- cut(old_comparison$ctr_copy_numbers, 
                                 summary(old_comparison$ctr_copy_numbers)[-4])
  old_comparison$de_status <- ifelse(old_comparison$abs_fold_change == 1,
                                     0, ifelse(
                                       old_comparison$abs_fold_change >= 1,
                                       1, -1))
  old_direction <- ifelse(old_comparison$b > 0, 1, -1)
  old_decision <- ifelse(old_comparison$de_status[old_filtered] == 0,
                         0, ifelse(
                           old_comparison$de_status[old_filtered] == 
                             old_direction[old_filtered],
                           1, 0))
  old_comparison$tpr <- old_comparison$fpr <- old_comparison$fdr <- NA
  old_comparison$tpr[old_filtered] <- cumsum(old_decision) / 
    sum(old_decision, na.rm=T)
  old_comparison$fpr[old_filtered] <- cumsum(
    -1*(old_decision-1)) /
    (-1*(sum(old_decision-1, na.rm=T)))
  old_comparison$fdr[old_filtered] <- cumsum(-1*(old_decision-1)) / 
    (seq(length(old_decision)))
  
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
  alr_comparison <- alr_comparison[order(alr_comparison$qval,
                                         alr_comparison$pval,
                                         alr_comparison$b), ]
  alr_filtered <- which(!is.na(alr_comparison$pval))
  alr_comparison$quartile <- cut(alr_comparison$ctr_copy_numbers, 
                                 summary(alr_comparison$ctr_copy_numbers)[-4])
  alr_comparison$de_status <- ifelse(abs(alr_comparison$alr_diff) < 1e-5,
                                     0, ifelse(
                                       alr_comparison$alr_diff > 0,
                                       1, -1))
  alr_direction <- ifelse(alr_comparison$b > 0, 1, -1)
  alr_decision <- ifelse(alr_comparison$de_status[alr_filtered] == 0,
                         0, ifelse(
                           alr_comparison$de_status[alr_filtered] == 
                             alr_direction[alr_filtered],
                           1, 0))
  alr_comparison$tpr <- alr_comparison$fpr <- alr_comparison$fdr <- NA
  alr_comparison$tpr[alr_filtered] <- cumsum(alr_decision) / 
    sum(alr_decision, na.rm=T)
  alr_comparison$fpr[alr_filtered] <- cumsum(
    -1*(alr_decision-1)) /
    (-1*(sum(alr_decision-1, na.rm=T)))
  alr_comparison$fdr[alr_filtered] <- cumsum(-1*(alr_decision-1)) / 
    (seq(length(alr_decision)))
  
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
  pdf(file.path(sleuth_dir, paste0(prefix, suffix)))
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
  pdf(file.path(sleuth_dir, paste0(prefix, suffix)))
  print(g2)
  invisible(dev.off())
  
  list(sleuth = old_comparison, alr = alr_comparison)
}

compare_sim_to_truth <- function(final_results, sleuth_dir = ".",
                                 de_probs, dir_probs,
                                 num_reps = c(10, 10), prefixes = NULL,
                                 gene_mode = FALSE) {
  group_ids <- c(rep(1, num_reps[1]), rep(2, num_reps[2]))
  if (is.null(prefixes)) {
    sim2truth <- lapply(seq_along(final_results$results), sim_to_truth,
                        sleuth_dir = sleuth_dir, gene_mode = gene_mode,
                        group_ids = group_ids, final_results = final_results)
  } else {
    sim2truth <- lapply(seq_along(final_results$results), function(x) {
      sim_to_truth(index = x, sleuth_dir = sleuth_dir, gene_mode = gene_mode,
                   group_ids = group_ids, final_results = final_results,
                   prefixes = prefixes[[x]])
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

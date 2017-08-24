# This function is the real meat of the simulation code
abs_simulation <- function(tpms, counts, s2c, eff_lengths,
                           sample_index = 1,
                           host = "dec2016.archive.ensembl.org",
                           species = "hsapiens", outdir = ".",
                           num_reps = c(10, 10),
                           gc_bias = rep(c(0, 0, 1, 1, 1, 2, 2, 3, 3, 3), 2),
                           de_prob = 0.01,
                           dir_prob = 0.5,
                           seed = 1,
                           polyester_sim = FALSE) {
  message("generating mean absolute copy numbers and relative TPMs")
  results <- generate_abs_changes(tpms = tpms,
                                  de_prob = de_prob,
                                  dir_prob = dir_prob,
                                  seed = seed,
                                  num_reps = num_reps)
  new_tpms <- results$transcript_abundances
  
  eff_lengths <- eff_lengths[match(rownames(new_tpms),
                                   eff_lengths$target_id), ]
  
  mean_lib_size <- 20*10^6
  
  message("converting relative TPMs to expected reads per transcript")
  expected_reads <- tpms_to_expected_reads(new_tpms,
                                           eff_lengths[,2],
                                           mean_lib_size)
  rel_fc <- expected_reads$fold_changes
  polyester_fc <- rel_fc
  polyester_fc[is.na(polyester_fc)] <- 1
  polyester_fc <- matrix(c(rep(1, length(polyester_fc)), polyester_fc),
                         nrow = length(polyester_fc))
  consistent <- calculate_consistency(results$abs_fold_changes,
                                      rel_fc)
  sizes <- NULL
  if (polyester_sim) {
    #require(polyester, lib.loc = "~/R/x86_64-pc-linux-gnu-library/3.4")
    lib_sizes <- rnorm(sum(num_reps), 1, sd = 0.01)
    reads_per_transcript <- expected_reads$expected_reads[, 1]#, drop = FALSE]
    sizes <- calculate_sizes(counts, s2c, reads_per_transcript,
                             polyester_fc[, 2])
    sizes[sizes==0] <- 1e-22
    #sizes = reads_per_transcript/3
    #sizes[sizes==0] <- 1e-22
    # how to empirically get GC bias
    #gc_bias <- ??
    print(length(reads_per_transcript))
    message("simulating an RNA-Seq experiment using 'polyester' package")
    polyester::simulate_experiment(fasta = fasta_file, outdir = outdir,
                                    num_reps = num_reps,
                                    reads_per_transcript = reads_per_transcript,
                                    size = sizes,
                                    fold_changes = polyester_fc,
                                    paired = TRUE, strand_specific = TRUE,
                                    distr = "empirical",
                                    error_model = "illumina5",
                                    bias = "none",
                                    gc_bias = gc_bias,
                                    #frag_GC_bias = gc_bias,
                                    seed = seed + 100000,
                                    gzip = TRUE, shuffle = TRUE)
  }
  return(append(results,
                values = list(sizes = sizes,
                              expected_reads = expected_reads$expected_reads,
                              adjusted_consistent_changes = consistent,
                              adjusted_fold_changes = rel_fc)))
}

#' Run a copy number simulation
#'
#' This function runs a simulation that explicitly defines copy numbers,
#' which can be used to test whether compositional changes leads to consistent
#' or misleading results based on the analysis done. After simulating an
#' experiment using the copy numbers, those numbers are converted into
#' expected reads to be used for a polyester simulation.
#' 
#' @param fasta_file, a multiFASTA file with the transcripts to be used in
#'   the simulation (required for polyester)
#' @param sleuth_file, a R-Data file containing the sleuth object containing
#'   results from a real experiment. If 'sleuth_save' is \code{FALSE}, the 
#'   sleuth object will be loaded using \code{load}, and the name of the
#'   sleuth object is expected to be 'sleuth.obj'.
#' @param sample_index, which sample from the real dataset should be used
#'   as the starting point for the simulation?
#' @param host, the URL to be used to download annotations using biomaRt.
#'   the default is the archive URL for Ensembl V87.
#' @param species, the abbreviated latin name of the species (default 
#'   "hsapiens"); used by biomaRt to get the annotations
#' @param outdir, where should the simulated reads be written to?
#' @param num_reps, the number of samples in each condition. Note that this
#'   only currently supports two conditions, so this must be length 2.
#' @param denom, the name(s) of transcript(s) that will be used as the
#'   denominator for the ALR transformation
#' @param seed, the random seed to be used for reproducibility
#' @param num_runs, the number of simulations to run
#' @param gc_bias, integer vector of length \code{sum(num_reps)} of the
#'   GC bias to be used by polyester. Only numbers between 0 and 7.
#'   See ?polyester::simulate_experiment under "gcbias" in the 
#'   Details section for more information.
#' @param de_probs, vector of same length as \code{num_runs}, with
#'   numbers between 0 and 1 describing the probability of differential
#'   expression for each simulation
#' @param dir_probs, vector of samel length as \code{num_runs}, with
#'   numbers between 0 and 1 describing the probability of differential
#'   expression being increased, given a transcript that is changing.
#' @param polyester_sim, should polyester be run? (default to \code{FALSE}
#'   to save time when you are merely interested in the ground truth)
#' @param sleuth_save if \code{TRUE}, the sleuth object was saved using
#'   'sleuth_save' and will be loaded using 'sleuth_load'.
#' @num_cores the number of cores to be used to run parallel simulations.
#'   the default is to use just one.
#' 
#' @return list with two members:
#'   + results: a list of lists, one entry for each simulation. Each
#'   simulation's results has the following entries:
#'   -- all of the entries returned by \code{generate_abs_changes}
#'   -- sizes: the size parameter for each transcript
#'   -- expected_reads: an N x 2 matrix with the expected number of 
#'      fragments for each transcript in each condition
#'   -- adjusted_consistent_changes: consistency comparing copy numbers
#'      to the relative data after normalization using the DESeq procedure
#'   -- adjusted_fold_changes: the fold changes perceived after normalization
#'      using the DESeq procedure
#'   + alr_data: a list of lists, one entry for each simulation. Each
#'   simulation's alr_data list contains the results from
#'   \code{calculate_rel_consistency}
#' @export
run_abs_simulation <- function(fasta_file, sleuth_file, sample_index = 1,
                               host = "dec2016.archive.ensembl.org",
                               species = "hsapiens", outdir = ".",
                               num_reps = c(10,10),
                               denom = "ENST00000396859.5",
                               seed = 1, num_runs = 1,
                               gc_bias = rep(c(0, 0, 1, 1, 1, 2, 2, 3, 3, 3),
                                             2),
                               de_probs = 0.1,
                               dir_probs = 0.5,
                               polyester_sim = FALSE,
                               sleuth_save = FALSE, num_cores = 1) {
  message("loading sleuth results object")
  if (sleuth_save)
    sleuth.obj <- sleuth::sleuth_load(sleuth_file)
  else
    load(sleuth_file)
  tpms <- sleuth:::spread_abundance_by(sleuth.obj$obs_raw, "tpm")
  counts <- sleuth:::spread_abundance_by(sleuth.obj$obs_raw, "est_counts")
  s2c <- sleuth.obj$sample_to_covariates
  ctr_samples <- which(s2c$condition == levels(s2c$condition)[1])
  s2c <- s2c[ctr_samples, ]
  counts <- counts[, ctr_samples]
  
  message("loading transcripts from FASTA file")
  transcripts <- Biostrings::readDNAStringSet(fasta_file)
  if(nchar(names(transcripts)[1]) > 15) {
    names(transcripts) <- sapply(names(transcripts), function(x) {
      strsplit(x, " ", fixed = T)[[1]][2]
    })
  } else {
    message("loading annotations from biomaRt")
    mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                             host = host,
                             dataset = paste(species, "gene_ensembl", sep = "_"))
    table <- biomaRt::getBM(mart = mart,
                            attributes = c("ensembl_transcript_id",
                                           "transcript_version"))
    table$target_id <- paste(table$ensembl_transcript_id,
                             table$transcript_version,
                             sep = ".")
    table$transcript_version <- NULL
    new_names <- table[match(names(transcripts), table$ensembl_transcript_id),
                       "target_id"]
    names(transcripts) <- new_names
  }
  stopifnot(all(grepl("^ENS", names(transcripts))))
  
  tpms <- tpms[names(transcripts), sample_index]
  if (sum(tpms) != 10^6) tpms <- tpms / sum(tpms) * 10^6
  
  eff_lengths <- unique(sleuth.obj$obs_raw[, c("target_id", "eff_len")])
  
  results <- vector(mode = "list", length = num_runs)
  rel_consistency <- vector(mode = "list", length = num_runs)
  results <- parallel::mclapply(seq(num_runs), function(i) {
    message(paste0("running run #", i))
    dir.create(file.path(outdir, paste0("run", i, "_fasta")), showWarnings = F)
    real_outdir <- file.path(outdir, paste0("run", i, "_fasta"))
    result <- abs_simulation(tpms, counts, s2c, eff_lengths,
                             sample_index, host, species, real_outdir,
                             num_reps, gc_bias, de_probs[i], dir_probs[i],
                             seed + (i-1)*5*10^5, polyester_sim)
    # polyester_files <- list.files(outdir, "^sample")
    # info_files <- list.files(outdir, "^sim")
    # old_files <- c(polyester_files, info_files)
    # new_files <- paste(paste0("run", i), old_files)
    # 
    # file_results <- file.rename(from = file.path(outdir, old_files),
    #                             to = file.path(outdir, new_files))
    result
  }, mc.cores = min(parallel::detectCores()-2, num_runs, num_cores))
  alr_data <- parallel::mclapply(seq(num_runs), function(i) {
    message(paste0("finding ALR consistency from run #", i))
    calculate_rel_consistency(results[[i]]$copy_numbers_per_cell,
                              results[[i]]$transcript_abundances,
                              denom = denom)
  }, mc.cores = min(parallel::detectCores()-2, num_runs, num_cores))
  return(list(results = results, alr_data = alr_data))
}

#' Summarize a group of copy number simulations
#'
#' This function summarizes for each simulation the consistency of the
#' fold changes comparing copy numbers to TPMs, and then the consistency
#' comparing the ALR transformation of the copy numbers versus that of 
#' the TPMs.
#'
#' @param abs_sim_runs_list, the output from \code{run_abs_simulation}
#' @param num_runs, the number of simulations
#' @param de_probs, vector of probabilities of differential expression
#' @param dir_probs, vector of probabilities of differential expression
#'   being increased, given a transcript that is changing
#'
#' @return a data frame with the following columns:
#'   + run_num: which simulation?
#'   + de_prob: probability of differential expression
#'   + dir_prob: probability of DE being up
#'   + rel*: the number of transcripts with each category of consistency
#'   between copy numbers and relative TPMs
#'   + rel_apparent_fc: the apparent fold change of relative TPMs of
#'   transcripts that are not changing
#'   + adj*: the number of transcripts with each category of consistency
#'   between copy numbers and DESeq-normalized counts
#'   + adj_apparent_fc: the apparent fold change of DESeq-normalized
#'   counts for transcripts that are not changing
#'   + alr: the number of transcripts that have consistent
#'   changes after ALR transformation of copy numbers vs relative TPMs
#'
#' @export
summarize_abs_sim_runs <- function(abs_sim_runs_list, num_runs,
                                   de_probs, dir_probs) {
  summaries <- data.frame()
  for (i in seq(num_runs)) {
    de_prob <- de_probs[i]
    dir_prob <- dir_probs[i]
    result <- abs_sim_runs_list$results[[i]]
    rel_consistency <- abs_sim_runs_list$alr_data[[i]]$alr_consistency
    rel_summary <- summary(result$consistent_changes)
    rel_z2c <- which(result$consistent_changes=="zero2change")
    rel_apparent_fc <- unique(round(result$rel_fold_changes[rel_z2c], 10))
    adj_summary <- summary(result$adjusted_consistent_changes)
    adj_z2c <- which(result$adjusted_consistent_changes=="zero2change")
    adj_apparent_fc <- unique(round(result$rel_fold_changes[rel_z2c], 10))
    alr_summary <- summary(rel_consistency)[3]
    rel_df <- data.frame(run_num = i, de_prob = de_probs[i],
                         dir_prob = dir_probs[i],
                         rel = rbind(rel_summary), rel_fc = rel_apparent_fc)
    adj_df <- data.frame(run_num = i, de_prob = de_probs[i],
                         dir_prob = dir_probs[i],
                         adj = rbind(adj_summary), adj_fc = adj_apparent_fc)
    alr_df <- data.frame(run_num = i, de_prob = de_probs[i],
                         dir_prob = dir_probs[i],
                         alr = alr_summary)
    summary_df <- merge(rel_df, adj_df,
                        by = c("run_num", "de_prob", "dir_prob"))
    summary_df <- merge(summary_df, alr_df,
                        by = c("run_num", "de_prob", "dir_prob"))
    summaries <- rbind(summaries, summary_df)
  }
  summaries
}
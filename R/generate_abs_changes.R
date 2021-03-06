#' Generate Absolute Changes
#'
#' This function simulates an experiment by making
#' changes to absolute copy numbers, using different
#' probabilities of differential expression, of 
#' direction of change, and of what level of differential
#' expression. It returns the copy numbers, relative TPMs,
#' the real fold changes and perceived fold changes of the
#' relative TPMs, as well as the consistency of these two
#' fold changes. NOTE: this only supports two condition
#' experiments at this time.
#'
#' @param tpms a vector of length equal to the number of targets.
#'   preferably named so that the results are tied to target names.
#' @param de_prob a single value greater than 0 and less than 1 to
#'   denote the probability of an individual target having differential
#'   expression. This is referring to the expected percentage of
#'   all features that will be differentially expressed (DE). If min_tpm
#'   is set, then this probability will be adjusted to make sure that
#'   the number of DE features among the filtered features matches the
#'   expected proportion of DE for all features.
#' @param dir_prob a single value greater than 0 and less than 1 to
#'   denote the probability of the differential expression being increased.
#' @param de_levels a numeric vector of the different possible levels of
#'   differential expression (default has three: a 'small' change of 25\%,
#'   a moderate change of 100\% / 2-fold, and a 'large' change of 300\% / 4-fold
#' @param seed what the seed is for the random number generator, so that
#'   the results are reproducible
#' @param num_reps an integer vector describing the number of replicates
#'   each condition will have.
#' @param min_tpm the minimum transcripts per million that determines
#'   which transcripts are expressed highly enough to be considered
#'   as potentially differentially expressed. This helps avoid simulating
#'   differential expression with low abundance transcripts. Note that this
#'   does not necessarily correspond to low-abundance transcripts if considering
#'   estimated counts. This can be set to \code{NULL}, \code{FALSE}, or 0 to
#'   turn filtering off and consider all features.
#'
#' @return a list with the following members:
#'   \itemize{
#'     \item abs_fold_changes: a vector of the absolute fold changes for each target
#'     \item rel_fold_changes: a vector of the perceived relative fold changes
#'       for each target
#'     \item consistent_changes: the consistency of fold changes as determined by
#'       \code{\link{calculate_consistency}}
#'     \item copy_numbers_per_cell: an N x 2 matrix of the expected copy numbers
#'       for each target in each condition
#'     \item transcript_abundances: an N x 2 matrix of the expected relative TPMs
#'       for each target in each condition
#'   }
#'
#' @importFrom stats rbinom
#' @importFrom truncnorm rtruncnorm
#' @export
generate_abs_changes <- function(tpms = NULL,
                                 de_prob = 0.1,
                                 dir_prob = 0.5,
                                 de_levels = c(1.25, 2, 4),
                                 de_type = "discrete",
                                 seed = 1,
                                 num_reps = c(10, 10),
                                 min_tpm = 1) {
  set.seed(seed)
  stopifnot(is.numeric(tpms))
  if(length(num_reps)!=2)
    stop("this currently only supports an experiment with two conditions")
  if (sum(tpms) != 10^6) tpms <- tpms/sum(tpms) * 10^6

  de_type <- match.arg(de_type, c("discrete", "normal"))
  if (de_type == "discrete" & is.null(de_levels)) {
    stop("if you are doing a discrete fold-change simulation, you must ",
         "specify the levels of differential expression")
  }

  if (is.null(min_tpm) || !min_tpm) min_tpm <- 0

  # The conceptual shift from relative TPMs to absolute copy numbers
  ctr_copy_numbers <- tpms
  tpm_filter <- which(ctr_copy_numbers >= min_tpm)
  eligible_trans <- ctr_copy_numbers[tpm_filter]
  # garbage collection to keep memory footprint small
  rm(tpms)
  gc()

  ## these are now treated as transcript copy numbers / cell
  num_trans <- length(eligible_trans)
  num_levels <- length(de_levels)

  adjusted_de_prob <- de_prob * length(ctr_copy_numbers) / num_trans
  if (adjusted_de_prob > 1) {
    warning(paste0("The specified 'de_prob' and 'min_tpm' leads to there being ",
                   "more expected differentially expressed features than the total ",
                   "number of filtered features. Setting adjusted DE prob to 1."))
    adjusted_de_prob <- 1
  }

  if (num_levels < 3 & de_type == "normal") {
    stop("if you are using a truncated normal, you must specify ",
         "at least three values for the 'de_levels' variable to set ",
         "the basement, the mean, and the standard deviation of the ",
         "distribution")
  }
  
  ## de_decision is a bernoulli trial to determine whether each transcript is
  ## differentially expressed or not
  de_decisions <- rbinom(num_trans, 1, adjusted_de_prob)
  de_hits <- which(de_decisions==1)
  num_de <- length(de_hits)
  ## dir_decision is a bernoulli trial to determine which direction
  ## differential expression will occur
  dir_decisions <- rep(NA, length(de_decisions))
  dir_decisions[de_hits] <- rbinom(num_de, 1, dir_prob)

  ## determine fold changes based on the 'coin flips' for
  ## differential expression, then the direction, then the level
  eligible_fcs <- sapply(seq_along(dir_decisions), function(x) {
    dir_decision <- dir_decisions[x]
    if(is.na(dir_decision)) return(1)
    if (de_type == "discrete") {
      index <- as.integer(cut(runif(1), seq(0, 1, 1/num_levels)))
      level <- de_levels[index]
    } else {
      level <- truncnorm::rtruncnorm(1, a = de_levels[1],
                                     mean = de_levels[2],
                                     sd = de_levels[3])
    }
    ifelse(dir_decision == 1, level, 1/level)
  })

  ## the "experimental condition" copy numbers are the
  ## control copy numbers * fold changes
  fold_changes <- rep(1, length(ctr_copy_numbers))
  fold_changes[tpm_filter] <- eligible_fcs
  exp_copy_numbers <- ctr_copy_numbers * fold_changes

  ## the absolute data is combining the control and experimental numbers
  ## these represent the "mean" values for each transcript in each condition
  abs_samples <- as.matrix(cbind(ctr_copy_numbers, exp_copy_numbers))
  ## by renormalizing the copy numbers to TPMs, it becomes
  ## relative/compositional data again
  relative_samples <- sweep(abs_samples, 2, colSums(abs_samples), "/") * 10^6
  colnames(relative_samples) <- c("ctr_tpms", "exp_tpms")
  ## calculate fold changes in the classic way by comparing the changes
  ## seen between the relative TPMs
  relative_fc <- relative_samples[,2] / relative_samples[,1]
  
  consistent <- calculate_consistency(fold_changes, relative_fc)
  return(list(abs_fold_changes = fold_changes,
              rel_fold_changes = relative_fc,
              consistent_changes = consistent,
              copy_numbers_per_cell = abs_samples,
              transcript_abundances = relative_samples))
}

#' importFrom utils data
add_spikeins <- function(results, spikein_mix = "Mix1", spikein_percent = 0.02) {
  stopifnot(is.numeric(spikein_percent) && length(spikein_percent) == 1)
  stopifnot(spikein_percent > 0 && spikein_percent < 1)
  stopifnot(is.character(spikein_mix))
  stopifnot(all(spikein_mix %in% c("Mix1", "Mix2")))

  curr_copy_nums <- results$copy_numbers_per_cell
  actual_percent <- (1 / (1 - spikein_percent)) - 1

  data(ERCC92_data, package = "absSimSeq")

  if (length(spikein_mix) > 2) {
    stop("'spikein_mix' has more than two elements. 'absSimSeq' currently only supports two condition experiments")
  } else if (length(spikein_mix) == 2) {
    spike_cols <- paste0(spikein_mix, "_molar_conc")
  } else {
    spike_cols <- rep(paste0(spikein_mix, "_molar_conc"), 2)
  }

  spikein_copy_nums <- as.matrix(ERCC92_data[, spike_cols])
  colnames(spikein_copy_nums) <- colnames(curr_copy_nums)

  ratio <- actual_percent * sum(curr_copy_nums[,1]) / sum(spikein_copy_nums[,1])
  spikein_copy_nums <- sweep(spikein_copy_nums, 2, ratio, "*")

  new_copy_nums <- rbind(curr_copy_nums, spikein_copy_nums)
  new_fold_changes <- new_copy_nums[,2] / new_copy_nums[,1]
  new_fold_changes[is.na(new_fold_changes)] <- 1

  new_tpms <- sweep(new_copy_nums, 2, colSums(new_copy_nums), "/") * 10^6
  new_rel_fcs <- new_tpms[,2] / new_tpms[,1]
  new_rel_fcs[is.na(new_rel_fcs)] <- 1

  consistent <- calculate_consistency(new_fold_changes, new_rel_fcs)

  return(list(abs_fold_changes = new_fold_changes,
              rel_fold_changes = new_rel_fcs,
              consistent_changes = consistent,
              copy_numbers_per_cell = new_copy_nums,
              transcript_abundances = new_tpms))
}


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
#'   expression
#' @param dir_prob a single value greater than 0 and less than 1 to
#'   denote the probability of the differential expression being increased.
#' @param de_levels a numeric vector of the different possible levels of
#'   differential expression (default has three: a "small" change of 25%,
#'   a moderate change of 100% / 2-fold, and a "large" change of 300% / 4-fold
#' @param seed what the seed is for the random number generator, so that
#'   the results are reproducible
#' @param num_reps an integer vector describing the number of replicates 
#'   each condition will have.
#' 
#' @return a list with the following members:
#'   + abs_fold_changes: a vector of the absolute fold changes for each target
#'   + rel_fold_changes: a vector of the perceived relative fold changes
#'     for each target
#'   + consistent_changes: the consistency of fold changes as determined by
#'     \link{calculate_consistency}
#'   + copy_numbers_per_cell: an N x 2 matrix of the expected copy numbers
#'     for each target in each condition
#'   + transcript_abundances: an N x 2 matrix of the expected relative TPMs
#'     for each target in each condition
#'
#' @export
generate_abs_changes <- function(tpms = NULL,
                                 de_prob = 0.1,
                                 dir_prob = 0.5,
                                 de_levels = c(1.25, 2, 4),
                                 seed = 1,
                                 num_reps = c(10, 10)) {
  set.seed(seed)
  stopifnot(is.numeric(tpms))
  if(length(num_reps)!=2)
    stop("this currently only supports an experiment with two conditions")
  if (sum(tpms) != 10^6) tpms <- tpms/sum(tpms) * 10^6
  
  # The conceptual shift from relative TPMs to absolute copy numbers
  ctr_copy_numbers <- tpms
  # garbage collection to keep memory footprint small
  rm(tpms)
  ## these are now treated as transcript copy numbers / cell
  num_trans <- length(ctr_copy_numbers)
  num_levels <- length(de_levels)
  
  ## de_decision is a bernoulli trial to determine whether each transcript is
  ## differentially expressed or not
  de_decisions <- rbinom(num_trans, 1, de_prob)
  de_hits <- which(de_decisions==1)
  num_de <- length(de_hits)
  ## dir_decision is a bernoulli trial to determine which direction
  ## differential expression will occur
  dir_decisions <- rep(NA, length(de_decisions))
  dir_decisions[de_hits] <- rbinom(num_de, 1, dir_prob)

  ## determine fold changes based on the 'coin flips' for
  ## differential expression, then the direction, then the level
  fold_changes <- sapply(seq_along(dir_decisions), function(x) {
    dir_decision <- dir_decisions[x]
    if(is.na(dir_decision)) return(1)
    index <- as.integer(cut(runif(1), seq(0, 1, 1/num_levels)))
    level <- de_levels[index]
    if(dir_decision==1) level else 1/level
  })
  
  ## the "experimental condition" copy numbers are the
  ## control copy numbers * fold changes
  exp_copy_numbers <- ctr_copy_numbers * fold_changes
  
  ## the absolute data is combining the control and experimental numbers
  ## these represent the "mean" values for each transcript in each condition
  abs_samples <- cbind(ctr_copy_numbers, exp_copy_numbers)
  ## by renormalizing the copy numbers to TPMs, it becomes
  ## relative/compositional data again
  relative_samples <- apply(abs_samples, 2, function(x) {
    x / sum(x) * 10^6
  })
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
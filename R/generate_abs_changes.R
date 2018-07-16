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
#'   differential expression (default has three: a 'small' change of 25\%,
#'   a moderate change of 100\% / 2-fold, and a 'large' change of 300\% / 4-fold
#' @param seed what the seed is for the random number generator, so that
#'   the results are reproducible
#' @param num_reps an integer vector describing the number of replicates
#'   each condition will have.
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

  if (num_levels < 3 & de_type == "normal") {
    stop("if you are using a truncated normal, you must specify ",
         "at least three values for the 'de_levels' variable to set ",
         "the basement, the mean, and the standard deviation of the ",
         "distribution")
  }
  
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

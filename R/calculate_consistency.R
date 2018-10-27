#' Calculate Consistency of Relative vs Absolute Changes
#'
#' If the assumption holds true that the overall
#' composition of (and the total amount of) RNA in each 
#' sample is the same, then the absolute and relative fold
#' changes should agree. If the assumption is violated,
#' then they will not agree. This function compares the 
#' absolute fold change and the perceived fold change
#' of the proportions of each target in a collection. It
#' reports a summary of how many targets were consistent
#' and how many were inconsistent. It also summarizes the
#' kind of consistency or inconsistency.
#'
#' @param fold_changes a vector of absolute fold changes
#' @param relative_fc a vector of perceived relative fold changes
#' @param log_scale if \code{TRUE} (default is \code{FALSE}),
#'   the fold_changes are on the log scale (i.e. 0 is no change);
#'   otherwise the fold_changes are on the absolute scale (i.e. 1
#'   is no change).
#' 
#' @return a factor vector listing the consistency for each target.
#'   the possible levels are "zeroCNs" (zero copy numbers),
#'   "consistent_noChange", "consistent_change", "zero2change"
#'   (an absolute 'no change' is perceived as a relative change),
#'   "change2zero" (an absolute change is perceived as not changing),
#'   "down2up" (an absolute decrease is perceived as a relative increase),
#'   or "up2down" (an absolute increase is perceived as a relative decrease)
#'
#' @export
calculate_consistency <- function(fold_changes = NULL, relative_fc = NULL,
                                  log_scale = FALSE) {
  if (log_scale) {
    result <- ifelse(is.na(relative_fc),
                     "zeroCNs",
                     ifelse((relative_fc > 0 &
                               fold_changes > 0) |
                              (relative_fc < 0 &
                                 fold_changes < 0),
                            "consistent_change",
                            ifelse(abs(relative_fc) < 10^-5 &
                                     fold_changes == 0,
                                   "consistent_noChange",
                                   ifelse(relative_fc > 0 &
                                            fold_changes < 0,
                                          "down2up",
                                          ifelse(relative_fc < 0 &
                                                   fold_changes > 0,
                                                 "up2down",
                                                 ifelse(relative_fc == 0 &
                                                          fold_changes != 0,
                                                        "change2zero",
                                                        "zero2change")
                                          )
                                   )
                            )
                     )
    ) 
  } else {
    result <- ifelse(is.na(relative_fc),
                     "zeroCNs",
                     ifelse((log(relative_fc) > 0 &
                               log(fold_changes) > 0) |
                              (log(relative_fc) < 0 &
                                 log(fold_changes) < 0),
                            "consistent_change",
                            ifelse(abs(log(relative_fc)) < 10^-5 &
                                     log(fold_changes) == 0,
                                   "consistent_noChange",
                                   ifelse(log(relative_fc) > 0 &
                                            log(fold_changes) < 0,
                                          "down2up",
                                          ifelse(log(relative_fc) < 0 &
                                                   log(fold_changes) > 0,
                                                 "up2down",
                                                 ifelse(log(relative_fc) == 0 &
                                                          log(fold_changes) != 0,
                                                        "change2zero",
                                                        "zero2change")
                                          )
                                   )
                            )
                     )
    ) 
  }
  factor(x = result, levels = c("zeroCNs", "consistent_noChange",
                               "consistent_change", "zero2change",
                               "change2zero", "down2up",
                               "up2down"))
}

#' Calculate Consistency of Logratio Data
#'
#' This function uses \link{calculate_consistency} to
#' determine if there is consistency in the perceived
#' fold changes after logratio transformation of the
#' absolute copy numbers and the relative TPMs. In
#' theory, all targets should show consistent behavior.
#'
#' @param abs_cns, an N x 2 matrix or data frame with expected
#'   absolute copy numbers in each of two conditions for N targets
#' @param rel_tpms, an N x 2 matrix or data frame with expected
#'   relative TPMs of N targets in each of two conditions
#' @param denom, the name(s) of the target(s) to be used as the
#'   denominator for the additive logratio transformation
#'
#' @return a list with three members:
#'   \itemize{
#'     \item alr_means: the ALR transformed values of the copy numbers
#'     \item alr_diff: the difference of logratios between the two conditions
#'       of the ALR transformed copy numbers
#'     \item alr_consistency: the result from \code{calculate_consistency}
#'   }
#' @importFrom sleuthALR alr_transformation
#' @export
calculate_rel_consistency <- function(abs_cns, rel_tpms, denom = NULL) {
  abs_cns <- as.matrix(abs_cns)
  if(is.null(denom)) {
    abs_fc <- abs_cns[,2] / abs_cns[,1]
    no_change_rows <- which(abs_fc == 1)
    denom <- rownames(abs_cns)[no_change_rows[1]]
  }
  rel_tpms <- as.matrix(rel_tpms)
  alr_cns <- sleuthALR::alr_transformation(abs_cns, denom)
  alr_cns_diff <- alr_cns[, 2] - alr_cns[, 1]
  alr_tpms <- sleuthALR::alr_transformation(rel_tpms, denom)
  alr_tpms_diff <- alr_tpms[, 2] - alr_tpms[, 1]
  list(alr_means = alr_cns, alr_diff = alr_cns_diff,
       alr_consistency = calculate_consistency(alr_cns_diff, alr_tpms_diff,
                                               log_scale = T),
       denom = denom)
}

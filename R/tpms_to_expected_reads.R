#' Convert TPMs to Expected Fragments
#' 
#' This function takes as input a matrix of transcript abundances
#' (in TPMs) and converts them to expected fragments for each
#' transcript in each condition. The expected reads is proportional
#' to the effective length for the transcript and the expected size
#' of the library.
#'
#' @param trans_abund, a N x 2 matrix or data frame of transcript
#'   abundances (in TPMs) for each transcript in each of two conditions.
#' @param eff_lens, a vector of effective lengths for each transcript.
#'   the transcripts must be in the same order as what is found in
#'   \code{trans_abund}.
#' @param mean_library_size, a single value of the expected size of the
#'   sequencing library (e.g. 60 million fragments)
#' 
#' @return a list with two entries:
#'   + expected_reads: an N x 2 matrix of the expected reads for each
#'   transcript in each condition
#'   + fold_changes: the perceived fold changes if one were to use the
#'   expected reads (unnormalized by library size or transcript effective
#'   length).
#' 
#' @export
tpms_to_expected_reads <- function(trans_abund = NULL,
                                   eff_lens = NULL,
                                   mean_library_size = NULL) {
  # translate tpms to read abundances --> 
  # then multiply by library size to get expected reads per transcript
  expected_counts <- sweep(trans_abund, 1, eff_lens, "*")
  read_abundances <- apply(expected_counts, 2, function(x) x / sum(x))
  expected_reads <- read_abundances * mean_library_size
  fold_changes <- expected_reads[, 2] / expected_reads[, 1]
  fold_changes[is.na(fold_changes)] <- 1
  size_factors <- DESeq2::estimateSizeFactorsForMatrix(counts = expected_reads)
  adjusted_reads <- sweep(expected_reads, 2, size_factors, FUN = "*")
  fold_changes <- adjusted_reads[, 2] / adjusted_reads[, 1]
  list(expected_reads = expected_reads,
       fold_changes = fold_changes)
}

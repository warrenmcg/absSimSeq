# This function takes a vector of IDs and their mean abundances (in TPMs)
# it returns only one ID, and prioritizes IDs with the highest expression
# and otherwise returns the first one with the lowest variation and
# the highest expression 
reduceToOneCandidate <- function(ids, means) {
  index <- NA
  if (length(ids)==1)
    index <- 1 
  else {
    id_means <- means(ids)
    max_means <- which(id_means == max(id_means))
    if(length(max_means) > 1) {
      max_mean_ids <- names(id_means[max_means])
      index <- which(names(ids) == max_mean_ids[1])
    } else {
      index <- which(names(ids) == names(max_means))
    }
  }
  ids[index]
}

# This function takes sleuth results, and attempts to identify a few categories
# of possible "optimal" reference genes to be used for the ALR transformation:
#
# candidate #1: the transcript with the least coefficient of variation
#   across all transcripts
# candidate #2: the transcript with the least coefficient of variation
#   across all transcripts in the highest 25% of mean expression
# candidate #3: the transcript with the least coefficient of variation
#   across all transcripts with a mean abundance of 10 TPM or
#   5000 counts/reads per base (usually top ~5-10%)
# if there are ties, it tries to prioritize the transcript with the highest
# abundance, and otherwise takes the first one.
get_least_var_targets <- function(sleuth.obj = NULL, which_var = "tpm") {
  if (is.character(sleuth.obj)) {
    tryCatch(sleuth.obj <- sleuth::sleuth_load(sleuth.obj),
             error = function(e) load(sleuth.obj))
  } else {
    stopifnot(is(sleuth.obj, "sleuth"))
  }

  which_var <- match.arg(which_var, c("tpm", "est_counts", "scaled_reads_per_base"))

  table <- sleuth:::spread_abundance_by(sleuth.obj$obs_norm, which_var)
  table <- table[which(rownames(table) %in% sleuth.obj$filter_df$target_id), ]
  
  cov <- apply(table, 1, function(x) sd(x) / mean(x))
  means <- apply(table, 1, mean)

  small_id <- names(which(cov == min(cov, na.rm=T)))

  topquart <- summary(means)[5]
  topquart_id <- names(which(cov == min(cov[which(means>=topquart)])))

  threshold <- ifelse(which_var == "tpm", 10, 5000)
  highexp_id <- names(which(cov == min(cov[which(means>=threshold)])))

  small_id <- reduceToOneCandidate(small_id, means)
  topquart_id <- reduceToOneCandidate(topquart_id, means)
  highexp_id <- reduceToOneCandidate(highexp_id, means)

  ids <- c(small_id, topquart_id, highexp_id)
  names(ids) <- c("small_id", "topquartile_id", "highTPM_id")
  ids
}

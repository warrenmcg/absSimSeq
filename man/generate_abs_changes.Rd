% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_abs_changes.R
\name{generate_abs_changes}
\alias{generate_abs_changes}
\title{Generate Absolute Changes}
\usage{
generate_abs_changes(tpms = NULL, de_prob = 0.1, dir_prob = 0.5,
  de_levels = c(1.25, 2, 4), de_type = "discrete", seed = 1,
  num_reps = c(10, 10), min_tpm = 1)
}
\arguments{
\item{tpms}{a vector of length equal to the number of targets.
preferably named so that the results are tied to target names.}

\item{de_prob}{a single value greater than 0 and less than 1 to
denote the probability of an individual target having differential
expression. This is referring to the expected percentage of
all features that will be differentially expressed (DE). If min_tpm
is set, then this probability will be adjusted to make sure that
the number of DE features among the filtered features matches the
expected proportion of DE for all features.}

\item{dir_prob}{a single value greater than 0 and less than 1 to
denote the probability of the differential expression being increased.}

\item{de_levels}{a numeric vector of the different possible levels of
differential expression (default has three: a 'small' change of 25\%,
a moderate change of 100\% / 2-fold, and a 'large' change of 300\% / 4-fold}

\item{seed}{what the seed is for the random number generator, so that
the results are reproducible}

\item{num_reps}{an integer vector describing the number of replicates
each condition will have.}

\item{min_tpm}{the minimum transcripts per million that determines
which transcripts are expressed highly enough to be considered
as potentially differentially expressed. This helps avoid simulating
differential expression with low abundance transcripts. Note that this
does not necessarily correspond to low-abundance transcripts if considering
estimated counts. This can be set to \code{NULL}, \code{FALSE}, or 0 to
turn filtering off and consider all features.}
}
\value{
a list with the following members:
  \itemize{
    \item abs_fold_changes: a vector of the absolute fold changes for each target
    \item rel_fold_changes: a vector of the perceived relative fold changes
      for each target
    \item consistent_changes: the consistency of fold changes as determined by
      \code{\link{calculate_consistency}}
    \item copy_numbers_per_cell: an N x 2 matrix of the expected copy numbers
      for each target in each condition
    \item transcript_abundances: an N x 2 matrix of the expected relative TPMs
      for each target in each condition
  }
}
\description{
This function simulates an experiment by making
changes to absolute copy numbers, using different
probabilities of differential expression, of 
direction of change, and of what level of differential
expression. It returns the copy numbers, relative TPMs,
the real fold changes and perceived fold changes of the
relative TPMs, as well as the consistency of these two
fold changes. NOTE: this only supports two condition
experiments at this time.
}

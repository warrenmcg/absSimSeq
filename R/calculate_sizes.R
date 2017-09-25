# This function calculates size factors for use with the polyester simulation
# using real data.
# In a negative binomial distribution, if the size is r, then the variance of
# the NB random variable is V = mean + mean^2 / r
# It is estimated using DESeq2's \code{estimateDispersions} function.
# Note: Dispersion = 1 / r
calculate_sizes <- function(counts = NULL, s2c, reads_pt = NULL,
                            fc = NULL, single_value = TRUE) {
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(counts),
                                        colData = s2c,
                                        design = formula("~1"))
  DESeq2::sizeFactors(dds) <- 1
  dds <- DESeq2::estimateDispersions(dds)
  if (single_value) {
    # size (i.e. r) = 1 / dispersion
    sizes <- 1 / DESeq2::dispersions(dds)
  } else {
    func <- dds@dispersionFunction
    basemeans <- matrix(c(reads_pt, reads_pt), nrow=length(reads_pt))
    basemeans[,2] <- basemeans[,1] * fc
    # size (i.e. r) = 1 / dispersion
    sizes <- 1 / func(basemeans)
  }
  sizes
}

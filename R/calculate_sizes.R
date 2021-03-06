#' This function calculates size factors for use with the polyester simulation
#' using real data.
#' In a negative binomial distribution, if the size is r, then the variance of
#' the NB random variable is V = mean + mean^2 / r
#' It is estimated using DESeq2's \code{estimateDispersions} function.
#' Note: Dispersion = 1 / r
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateSizeFactors estimateDispersions
#'   dispersions dispersionFunction counts
#' @importFrom GenomicRanges mcols
#' @importFrom stats median
calculate_sizes <- function(counts = NULL, s2c, design, reads_pt = NULL,
                            fc = NULL, single_value = TRUE,
                            min_dispersion = 1e-6) {
  mode(counts) <- "integer"
  dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                        colData = s2c,
                                        design = design)
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- suppressMessages(DESeq2::estimateDispersions(dds))

  if (single_value) {
    dispersions <- DESeq2::dispersions(dds)
    basemeans <- data.frame(basemean = GenomicRanges::mcols(dds)$baseMean)
    rownames(basemeans) <- names(dispersions) <- rownames(DESeq2::counts(dds))
  } else {
    func <- DESeq2::dispersionFunction(dds)
    basemeans <- matrix(c(reads_pt, reads_pt), nrow=length(reads_pt))
    basemeans[, 2] <- basemeans[, 1] * fc
    dispersions <- func(basemeans)
    rm(func)
  }
  rm(dds)
  gc()

  disp_filt <- dispersions > min_dispersion
  disp_filt <- ifelse(is.na(disp_filt), FALSE, disp_filt)
  if (is.null(dim(dispersions))) {
    med_dispersion <- median(dispersions[disp_filt])
    dispersions[!disp_filt] <- med_dispersion
  } else {
    med_dispersions <- sapply(1:2, function(i) {
      median(dispersions[, i][disp_filt[, i]])
    })
    dispersions <- sapply(1:2, function(i) {
      med_disp <- med_dispersions[i]
      disps <- dispersions[, i]
      disps[!disp_filt[, i]] <- med_disp
      disps
    })
    dispersions <- t(dispersions)
    stopifnot(dim(dispersions) == dim(basemeans))
  }  
  # size (i.e. r) = 1 / dispersion
  sizes <- 1 / dispersions
  list(sizes = sizes, means = basemeans)
}

#' @importFrom stats median
calculate_spikein_sizes <- function(deseq_res, reads_pt, func = median) {
  stopifnot(!is.null(dim(reads_pt)))
  spikein_means <- rowMeans(reads_pt)[grepl("^ERCC", rownames(reads_pt))]
  gene_means <- deseq_res$means
  spike_sizes <- sapply(names(spikein_means), function(id) {
    mean <- spikein_means[id]
    matches <- which(signif(gene_means, 1) == signif(mean, 1))
    if (length(matches) == 0) {
      matches <- which(10^ceiling(log10(gene_means)) == 10^ceiling(log10(mean)))
      if (length(matches) == 0) {
        matches <- which(10^floor(log10(gene_means)) == 10^floor(log10(mean)))
        if (length(matches) == 0) {
          matches <- 1:length(gene_means)
        }
      }
    }
    matched_sizes <- deseq_res$sizes[matches]
    func(matched_sizes)
  })
  rm(func)
  spike_sizes
}


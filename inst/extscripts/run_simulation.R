fasta_file <- "~/warren/Annotations/HsEnsembl_v87/Genes/gencodeTranscripts.v25.fa"
sleuth_file <- "~/warren/re-analysis/sleuthAnalysis.RData"
denom <- "ENST00000396859.5" # GAPDH;
outdir <- "polyester_sim"
sample_index <- 3 # this is a control sample from the real dataset
num_runs <- 9
num_reps <- c(5,5)
gc_bias <- rep(c(0, 0, 1, 1, 1), 2)
de_probs <- c(0.01, 0.01, 0.01, 0.1, 0.1, 0.1, 0.25, 0.25, 0.25)
dir_probs <- c(0.5, 0.9, 0.1, 0.5, 0.9, 0.1, 0.5, 0.9, 0.1)

final_results <- run_abs_simulation(fasta_file, sleuth_file, sample_index,
                                    outdir = outdir, denom = denom,
                                    num_runs = num_runs, num_reps = num_reps,
                                    de_probs = de_probs, dir_probs = dir_probs,
                                    polyester_sim = T, gc_bias = gc_bias,
                                    num_cores = num_runs)
                                    #num_cores = num_runs)

summaries <- summarize_abs_sim_runs(final_results, num_runs, de_probs, dir_probs)
save(final_results, summaries, sample_index, de_probs, dir_probs, fasta_file,
     sleuth_file, outdir, gc_bias, num_runs, num_reps, denom,
     file = file.path(outdir, "polyester_groundTruth.RData"))

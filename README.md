# absSimSeq

**NOTE**: this tool is in alpha stage of development. If you choose to try it, any constructive feedback is
greatly appreciated! A vignette is coming soon!

This package simulates RNA-Seq experiments with an explicit definition of the 
absolute RNA counts.

This package extends the [`polyester`](https://github.com/alyssafrazee/polyester)
package by simulating an RNA-Seq using absolute RNA counts before
converting the result into relative proportions (in TPMs and expected
fragments) that are then used as input for `polyester`.

The basic workflow is the following:
1) Take an estimate of TPM values from a sleuth object.
2) Conceptually shift the unit from TPM to copy numbers per cell.
3) Simulate differential expression on the copy numbers
4) Convert the new copy numbers back to TPMs
5) Convert TPM values to expected reads per transcript
6) Simulate reads using the *polyester* ([bioconductor page](http://bioconductor.org/packages/release/bioc/html/polyester.html)) package.

## Introduction to Simple Command

To run the simulation based on copy numbers, run the following command:

```
final_results <- run_abs_simulation(fasta_file, sleuth_file, sample_index,
                                    outdir = outdir, denom = denom,
                                    num_runs = num_runs, num_reps = num_reps,
                                    de_probs = de_probs, dir_probs = dir_probs,
                                    polyester_sim = T, gc_bias = gc_bias,
                                    num_cores = num_runs)
```

After you run `kallisto` and then `sleuth` on the simulations, you can then 
run the following command to compare the results to the truth:

```
comparisons <- compare_sim_to_truth(final_results, sleuth_dir = outdir,
                                    de_probs = de_probs, dir_probs = dir_probs,
                                    num_reps = num_reps)
```

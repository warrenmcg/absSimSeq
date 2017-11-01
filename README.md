# absSimSeq

Simulate RNA-Seq experiments with an explicit definition of the 
absolute RNA counts.

This package extends the [`polyester`](https://github.com/alyssafrazee/polyester)
package by simulating an RNA-Seq using absolute RNA counts before
converting the result into relative proportions (in TPMs and expected
fragments) that are then used as input for `polyester`.

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

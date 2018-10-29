# absSimSeq

**NOTE**: A preprint is coming soon! If you choose to try this, any constructive feedback is
greatly appreciated!

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

## Usage

See the [vignette](https://github.com/warrenmcg/absSimSeq/blob/master/vignettes/absSimSeq.Rmd) for full details!

# Trans-eQTL
Trans-eQTL joint with Barbara Engelhardt

# Directory Structure

The directory structure here, and features of the `analysis` subdirectory (including the `Makefile`), are based on
[https://github.com/jdblischak/singleCellSeq](https://github.com/jdblischak/singleCellSeq).

Here's a brief summary of the directory structure.
- analysis: Rmd files for investigations; will generate figures in `figure/` subdirectory
- R: R scripts/functions used in analysis; no code that is actually run put here
- output: largish files of data, typically created by code (e.g. post-processed data, simulations)
- code: code used for preprocessing, generation of output files etc ; may be long-running
- data: datasets that generally will not change once deposited
- paper: the paper
- packrat: a directory that contains information about all the R package used.
See the R package `packrat` for more details.
- talks: any presentations
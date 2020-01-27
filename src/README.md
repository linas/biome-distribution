
Data processing
---------------
Some of the primary files:

* `load-files.scm` - Load genomic/proteomic datasets into the AtomSpace
* `cleanup-data.scm` - Scrub datasets for cleaner results.
* `bio-loop.scm` - Data-mine triangular and pentagonal interactions.
* `data-export.scm` - Export selected stats to CSV files for graphing.
* `gene-pairs.scm` - Provide a sparse-matrix API for data analysis.

Some of the secondary files:
* `performance-instrumentation.scm` - measure performance stats.
* `debug-wrapper.scm` - Used for debugging the MOZI code-base.
* `lmpd-genes.scm` - a list of 1482 genes.
* `mem-use.sh` - log RAM usage statistics.
* `grid.conf` and `path.conf` - handy-dandy cogserver configs

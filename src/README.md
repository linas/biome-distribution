
Data processing
---------------
Some of the primary files:

* `load-files.scm` - Load genomic/proteomic datasets into the AtomSpace
* `cleanup-data.scm` - Scrub datasets for cleaner results.
* `bio-loop.scm` - Data-mine triangular and pentagonal interactions.
* `bio-tetra.scm` - Data-mine tetrahedron interactions.
* `data-export.scm` - Export selected stats to CSV files for graphing.
* `gene-pairs.scm` - Provide a sparse-matrix API for data analysis.

The primary data-processing pipeline is to start up a guile shell,
and then load each of the above files in turn.  Review thier contents,
and perform the "obvious" operations that each of them suggest.
See also the `graphs` directory for assorted data-graphing scripts.

Some of the secondary files:
* `performance-instrumentation.scm` - measure performance stats.
* `debug-wrapper.scm` - Used for debugging the MOZI code-base.
* `lmpd-genes.scm` - a list of 1482 genes.
* `mem-use.sh` - log RAM usage statistics.
* `grid.conf` and `path.conf` - handy-dandy cogserver configs

Typical run
-----------
A tetrahedron run:

```
HUGETLB_MORECORE=yes LD_PRELOAD=/usr/lib/libhugetlbfs.so.0 guile
scheme@(guile-user)> (load "load-files.scm") ; about 100 seconds
(start-cogserver)
(load-all)
(load "cleanup-data.scm")
(delete-go-nodes)
(delete-self-interaction)
(delete-bad-chebi)
(count-gene-interactions)
(symmetrize-gene-interactions)
(count-gene-interactions)
(delete-simple-tv)     ; about 160 seconds
(load "bio-loop.scm")
(make-triangles)
; Above should create 1797281 (1.8M) triangles
(load "bio-tetra.scm")
(count-tetrahedra (cog-get-atoms 'GeneNode))
```

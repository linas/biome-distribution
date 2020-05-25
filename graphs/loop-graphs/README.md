
Interaction graph distribution
------------------------------
Obtain distribution of gene interactions. 

See the email chain titled “Genome distribution!?” and also
the opencog benchmark “query-loop” for details about distributions
of interactions.

### HOWTO:

```
HUGETLB_MORECORE=yes LD_PRELOAD=/usr/lib/libhugetlbfs.so.0 guile
scheme@(guile-user)> (load "load-files.scm")
(start-cogserver)
(load-all)             ; about 100 seconds
(load "cleanup-data.scm")
(delete-go-nodes)
(delete-self-interaction)
(delete-bad-chebi)
(count-gene-interactions)
(make-gene-pairs)      ; 20 seconds
(cog-delete-recursive (Predicate "interacts_with"))
(delete-simple-tv)     ; about 160 seconds
(load "bio-loop.scm")
(count-gene-triangles (cog-get-atoms 'GeneNode))   ; about 2800 seconds
; Above should create 1797281 (1.8M) triangles
```

At this time, after deletetion of self-interactors, have
```
Found 537839 gene interactions
Found 344188 symmetric (paired) gene interactions
Conclude: there are 193651 asymmetric interctions
Conclude: there are 365745 symmetrized interactions
```

For the pentagons,
```
(define pathways (pathways-of-genes (cog-get-atoms 'GeneNode)))
(count-pentagons pathways)
```

To make the csv's, do a `(load "data-export.scm")` and then start
cutting and pasting the code there.

### Files

* `gene-loops.gplot` -- Create Zipf-distribution graphs for triangles.
* `tri-edges.gplot` -- Create Zipf-distribution graphs for gene-pair edges in
                       appearing triangles

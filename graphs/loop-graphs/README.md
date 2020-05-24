
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
(make-gene-pairs)
(cog-delete-recursive (Predicate "interacts_with"))
(delete-simple-tv)     ; about 160 seconds
(load "bio-loop.scm")
(count-gene-triangles)   ; about 2800 seconds
; Above should create 1797281 (1.8M) triangles
```


Then see bottom of that file...
```
(count-triangles (cog-get-atoms 'GeneNode))
; or
(define pathways (pathways-of-genes (cog-get-atoms 'GeneNode)))
(count-pentagons pathways)
```

To make the csv's, do a `(load "gene-pairs.scm")` and then go to
`data-export.scm` and start cutting and pasting the code there.

### Files

* `gene-loops.gplot` -- Create Zipf-distribution graphs for triangles.
* `tri-edges.gplot` -- Create Zipf-distribution graphs for gene-pair edges in
                       appearing triangles

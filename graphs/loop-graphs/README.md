
Interaction graph distribution
------------------------------
Obtain distribution of gene interactions. 

See the email chain titled “Genome distribution!?” and also
the opencog benchmark “query-loop” for details about distributions
of interactions.

### HOWTO:

```
(load "load-files.scm")
(start-cogserver)
(load-all)
(load "bio-loop.scm")
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

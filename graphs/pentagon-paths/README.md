
Pentagonal pathway-protein-gene distribution
--------------------------------------------
Obtain distribution of genes interacting, that express proteins that
interact on a common pathway.

See the [../loop-graphs/README.md](../loop-graphs/README.md) file first.

### HOWTO:

```
(load "load-files.scm")
(start-cogserver)
(load-all)
(load "bio-loop.scm")
```
Then see bottom of that file...
```
(define pathways (pathways-of-genes (cog-get-atoms 'GeneNode)))
(count-pentagons pathways)
```

To make the csv's, do a `(load "gene-pairs.scm")` and then go to 
`data-export.scm` and start cutting and pasting the code there.

### Files

* `path-vertex.gplot` -- Create Zipf-distribution graphs for vertexes of pentagons.
* `path-vertex-asym.gplot` -- As above, but for the older, faulty dataset.
                              Sanity check.
* `path-edges.gplot` -- Create Zipf-distribution graphs for edges of pentagons.

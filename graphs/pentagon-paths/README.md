
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

To make the csv's, cut and paste the code at the bottom of that file.

### Files

* `path-distribution.gplot` -- Create Zipf-distribution graphs for pentagons.

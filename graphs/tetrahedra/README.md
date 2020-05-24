
Tetrahedron graph distribution
------------------------------
Obtain distribution of tetrahedral gene interactions. 

See the README in `loop-graphs` for earlier discussion.

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
(symmetrize-gene-interactions)
(count-gene-interactions)
(delete-simple-tv)     ; about 160 seconds
(load "bio-loop.scm")
(make-triangles)       ; about 2800 seconds
; Above should create 1797281 (1.8M) triangles
(load "bio-tetra.scm")
(count-tetrahedra (cog-get-atoms 'GeneNode))
```

To make the csv's, do a `(load "data-export.scm)` and then start
cutting and pasting the code there.

Current results:
9989 genes participated in tetrahedra - `(length graph-participants)`

### Files

* `gene-tetra.gplot` -- Create Zipf-distribution graphs for tetrahedra.
* `tetra-edges.gplot` -- Create Zipf-distribution graphs for gene-pair
                        edges in appearing tetrahedra.

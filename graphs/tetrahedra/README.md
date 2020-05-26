
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
; (count-gene-interactions)
(make-gene-pairs)  ; about 20 seconds
; Above should make 365745 pairs
(delete-all-but-gene-interactions) ; about 66 seconds
(delete-simple-tv)     ; about 15 seconds
(load "bio-loop.scm")
(make-gene-triangles)       ; about 4970 seconds
; Above should create 1797281 (1.8M) triangles
(load "bio-tetra.scm")
(count-tetrahedra (cog-get-atoms 'GeneNode))
```

To make the csv's, do a `(load "data-export.scm)` and then start
cutting and pasting the code there.

Current results:
9989 genes participated in tetrahedra - `(length graph-participants)`
Theese were counted 884761344 but this is a 16x over-count, so actually
just 55297584 counts = 55M wow.

487728 distinct directed edges in tetrahedra - `(length gene-pairs)`
       (edges with non-zero counts) 
243864 symmetric edges in tetrahedra (half of above)
1797281 triangles, as always, but .. we didn't count triangles. :-(

1327142016.0 total observations of edges


### Files

* `gene-tetra.gplot`  -- Create Zipf-distribution graphs for tetrahedra.
* `tetra-edges.gplot` -- Create Zipf-distribution graphs for gene-pair
                         edges in appearing tetrahedra.

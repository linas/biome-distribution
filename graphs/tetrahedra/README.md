
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
(load "data-export.scm")
(print-loop-report)
```

To make the csv's, cut and paste from "data-export.scm".

Current results:
9989 genes participated in tetrahedra - `(length graph-participants)`
These were counted 147460224 times but this is a 4x over-count,
so actually 36865056 corner observations. 
Each tetra has 4 corners so 36865056 / 9216264 = 4 yay!

Total of 243864 symmetric edges participated in tetrahedra 
(out of 365745 edges total in the dataset.)
These were observed 221190336 times but this is a 4x overcount ...
4x because each tetrahedron was observed 4x
So total of 55297584 edge observations.
Each tetra has 6 edges so 55297584 / 9216264 = 6 yayy!

Participating triangles are 1701579 out of 1797281 triangles total
Observed 147460224 times, a 4x over-count.
4x because each tetrahedron was observed 4x
So a total of 36865056 triangle observations.
Each tetra has 4 tris so 36865056 / 9216264 = 4 yay!

Total of 9216264 tetrahedra, observed 36865056 times or a 4x over-count.
which makes sense: each was counted once, per corner.

### Files

* `gene-tetra.gplot`  -- Create Zipf-distribution graphs for tetrahedra.
* `tetra-edges.gplot` -- Create Zipf-distribution graphs for gene-pair
                         edges in appearing tetrahedra.

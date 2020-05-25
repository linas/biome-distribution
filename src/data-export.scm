;
; data-export.scm
;
; Extract counts from the AtomSpace graph, dump to CSV files for graphing.
;
(use-modules (srfi srfi-1))
(use-modules (opencog) (opencog exec))

; Needeed for definition of GeneNode and MoleculeNode
(use-modules (opencog bioscience))

; -----------------------------------------------------------------

(define (dump-to-csv pair-list filename)
"
   Write the pair-list to the filename.
   pair-list is a list of (string . count) pairs.
   Before writing, it will be sorted, first.
"
	; Sort according to descending rank.
	(define sorted-counts (sort pair-list
		(lambda (a b) (> (cdr a) (cdr b)))))

	; Dump to file.
	(define f (open-file filename "w"))
	(define cnt 1)
	(format f "#\n# ~A\n#\n# Rank-ordered counts\n#\n" filename)
	(for-each
		(lambda (gu) (format f "~A	~A	~A\n" cnt (car gu) (cdr gu))
			(set! cnt (+ 1 cnt)))
		sorted-counts)
	(close f)
)

#! -----------------------------------------------------------------
; Some stuff to create a ranked graph of the results found above.
; This first section is for the triangles; for the pentagons, see
; further down.

; Grand total number of participating genes.
; For stand-alone edges, this should be 20123.
; For triangles, this is 13846.
; For tetrahedra, this is 9989.
(define total-genes
	(fold (lambda (gene cnt) (if (< 0 (cog-count gene)) (+ cnt 1) cnt)) 0 
		(cog-get-atoms 'GeneNode)))

; Grand total observation count on all genes. i.e. how many times that
; gene was counted as participating in some graph.
;
; When counting edges, this should be 4x the number of edges,
; or 4x 365745 so 1462980. Why 4x? There is 2x because each endpoint is
; explored, and another 2x because UnorderedLink permutes.
;
; When counting triangles, this is 1462980 which is a 12x over-count:
; 3x because once per corner, 2x because distal edge is swapped, and
; another 2x because UnorderedLink explores both permutations. So the
; actual count is 1462980 / 12 = 121915
;
; When counting tetrahedra, this is 884761344
(define total-gene-observation-count
	(fold (lambda (gene cnt) (+ cnt (cog-count gene))) 0 
		(cog-get-atoms 'GeneNode)))

; Genes that had postive counts.
; That is, they were counted as participating in a graph.
(define graph-participants
	(map (lambda (gene) (cons (cog-name gene) (cog-count gene)))
		(filter (lambda (gene) (< 0 (cog-count gene)))
			(cog-get-atoms 'GeneNode))))

; Triangle counts go into "gene-loops.csv", tetrahedra into "gene-tetra.csv"
(dump-to-csv graph-participants "gene-pairs.csv")
(dump-to-csv graph-participants "gene-loops.csv")
(dump-to-csv graph-participants "gene-tetra.csv")

; ------------------------------------------------------------
; Above was for single vertexes. Below is for edges.
;
; Gene pairs that appeared as edges in a graph (triangle, tetrahedron
; or other graph).
; Some of these pairs have AnyNodes from the matrix code,
; so filter those out... 
; Some pairs have zero-count. These are interacting genes in 
; the dataset, but are not graph participants.
;
; For triangles, the length of this is 308765 (out of a total of 365745
; gene-pairs in the dataset).
;
(define participating-gene-pairs
	(filter (lambda (evlnk)
			(and (< 0 (cog-count evlnk))
				(equal? 'GeneNode (cog-type (gadr evlnk)))
				(equal? 'GeneNode (cog-type (gddr evlnk)))))
		(cog-incoming-set (Predicate "gene-pair"))))

; Total count on all of the edges participating in a graph.
; For triangles, this was 16175529, which is 9x the number of unique
; triangles. Why 9x? A factor of 3x because each vertex is visted.
; That is, each triangle is counted 3 times. A second factor of 3x
; because each edge is counted once. (!) We didn't really need to count
; each edge, by symmetry they would have shown up, eventually.
(define total-edge-observation-count
	(fold (lambda (pare cnt) (+ cnt (cog-count pare))) 0 
		(cog-incoming-set (Predicate "gene-pair"))))

; Count-pairs for the gene-pairs
(define gene-pair-cnts
	(map (lambda (evelnk) (cons
		(string-concatenate
			(list (cog-name (gadr evelnk)) "-" (cog-name (gddr evelnk))))
		(cog-count evelnk)))
		participating-gene-pairs))

; Use "tri-edges.csv" for triangles, "tetra-edges.csv" for tetrahedra.
(dump-to-csv gene-pair-cnts "tri-edges.csv")
(dump-to-csv gene-pair-cnts "tetra-edges.csv")

; ------------------------------------------------------------
; Above was for edges.  This is for triangles.
; Avoid errors:
(for-each cog-delete-recursive (cog-get-atoms 'Variable))
(cog-report-counts)

; How many triangles?
(define total-triangles
	(cog-incoming-size (Predicate "gene-triangle")))

; How many triangles were NOT observed exactly 3 times?
; Should be zero for plain triangles.
; but something else for tetrahedra.
(length (filter (lambda (evlnk) (not (= 3 (cog-count evlnk))))
	(cog-incoming-set (Predicate "gene-triangle"))))

; ------------------------------------------------------------
; Examination of the pentagons.

; Genes that appeared in a pentagonal loop.
(define path-genes
	(map (lambda (gene) (cons (cog-name gene) (cog-count gene)))
		(filter (lambda (gene) (< 0 (cog-count gene)))
			(cog-get-atoms 'GeneNode))))

(dump-to-csv path-genes "path-genes.csv")

; Proteins that appeared in a pentagonal loop.
(define path-proteins
	(map (lambda (protein) (cons (cog-name protein) (cog-count protein)))
		(filter (lambda (protein) (< 0 (cog-count protein)))
			(cog-get-atoms 'MoleculeNode))))

(dump-to-csv path-proteins "path-proteins.csv")

; Pathways that appeared in a pentagonal loop.
(define path-loops
	(map (lambda (pathway) (cons (cog-name pathway) (cog-count pathway)))
		(filter (lambda (pathway) (< 0 (cog-count pathway)))
			(cog-get-atoms 'ConceptNode))))

(dump-to-csv path-loops "path-loops.csv")

; Path-protein edges that appeared in a pentagonal loop.
(define path-edges
	(map (lambda (memb) (cons
			(string-append (cog-name (cog-outgoing-atom memb 0)) "-x-"
				(cog-name (cog-outgoing-atom memb 1)))
		(cog-count memb)))
		(filter (lambda (memb) (< 0 (cog-count memb)))
			(cog-get-atoms 'MemberLink))))

(dump-to-csv path-edges "path-edges-sym.csv")

; Protein-expessed-by-gene edges
(define path-exprs
	(map (lambda (expr) (cons
			(string-append (cog-name (gadr expr)) "-x-"
				(cog-name (gddr expr)))
		(cog-count expr)))
		(filter (lambda (expr)
				(and (< 0 (cog-count expr))
					; Avoid inclusion of 'AnyNode
					(eq? (cog-type (gadr expr)) 'GeneNode)
					(eq? (cog-type (gddr expr)) 'MoleculeNode)))
			(cog-incoming-by-type (Predicate "expresses") 'EvaluationLink))))

(dump-to-csv path-exprs "path-exprs-sym.csv")

; Gene interactions
(define path-intrs
	(map (lambda (intr) (cons
			(string-append (cog-name (gadr intr)) "-x-"
				(cog-name (gddr intr)))
		(cog-count intr)))
		(filter (lambda (intr) (< 0 (cog-count intr)))
			(cog-incoming-by-type (Predicate "interacts_with") 'EvaluationLink))))

(dump-to-csv path-intrs "path-intrs-sym.csv")

; How many pathways? Lets count. I get 2129.
(fold (lambda (path cnt) (+ cnt (if (< 0 (cog-count path)) 1 0))) 0
	(cog-get-atoms 'ConceptNode))

; How many pentagons? Lets count paths. I get 491558.0
(fold (lambda (path cnt) (+ cnt (cog-count path))) 0
	(cog-get-atoms 'ConceptNode))

; How many pentagons, via proteins? I get 983116.0 twice as many, of course
(fold (lambda (prot cnt) (+ cnt (cog-count prot))) 0
	(cog-get-atoms 'MoleculeNode))

; Via genes? 983116.0 - twice as many: its mirror symmetry.
(fold (lambda (gene cnt) (+ cnt (cog-count gene))) 0
	(cog-get-atoms 'GeneNode))

; Proteins belonging to pathways? 983116.0
(fold (lambda (memb cnt) (+ cnt (cog-count memb))) 0
	(cog-get-atoms 'MemberLink))

; The other three? 2457790.0 = 2 * 983116.0 + 491558.0  OK
(fold (lambda (eval cnt) (+ cnt (cog-count eval))) 0
	(cog-get-atoms 'EvaluationLink))

; ------------------------------------------------------------------
; Gene-expression distribution.

(define gea (make-expression-pair-api))
(define ges (add-pair-stars gea))
(define gez (add-zero-filter ges #f))
(gez 'left-basis-size) ; 6694
(gez 'right-basis-size) ; 6735
(batch-all-pair-mi gez)
(print-matrix-summary-report gez)

; Distribution of support -- how many proteins each gene expresses
(define expr-supp
	(map (lambda (gene) (cons (cog-name gene)
		(fold (lambda (PR cnt)
			(if (< 0 (cog-count PR)) (+ 1 cnt) cnt)) 0
			(gez 'right-stars gene))))
		(gez 'left-basis)))

(dump-to-csv expr-supp "expr-supp.csv")

; Should be identical to path-genes, up above.
(define expr-cnt
	(map (lambda (gene) (cons (cog-name gene)
		(fold (lambda (PR cnt) (+ cnt (cog-count PR))) 0
			(gez 'right-stars gene))))
		(gez 'left-basis)))

(dump-to-csv expr-cnt "expr-cnt.csv")

!# ; ---------------------------------------------------------------

; ------------------------------------------------------------------

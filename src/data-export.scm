;
; data-export.scm
;
; Extract counts from the AtomSpace graph, dump to CSV files for graphing.
;
(use-modules (srfi srfi-1))
(use-modules (opencog) (opencog exec))

; Needeed for definition of GeneNode and MoleculeNode
(use-modules (opencog bioscience))

; =================================================================
; Print summary reports.

(define (print-loop-report)
	; Avoid errors:
	(for-each cog-delete-recursive (cog-get-atoms 'Variable))
	(cog-report-counts)

	; Grand total number of participating genes.
	; For stand-alone edges, this should be 20123.
	; For triangles, this is 13846.
	; For tetrahedra, this is 9989.
	(define total-genes
		(fold (lambda (gene cnt) (if (< 0 (cog-count gene)) (+ cnt 1) cnt)) 0
			(cog-get-atoms 'GeneNode)))
	(format #t "Found ~A genes with non-zero count, out of ~A total\n"
		total-genes (length (cog-get-atoms 'GeneNode)))

	; -------
	; Grand total observation count on all genes. i.e. how many times that
	; gene was counted as participating in some graph.
	;
	; When counting edges, this should be 4x the number of edges,
	; or 4x 365745 so 1462980. Why 4x? There is 2x because each endpoint is
	; explored, and another 2x because UnorderedLink permutes.
	;
	; When counting triangles, this is 16175529 which is a 9x over-count.
	; 3x because once per corner, and 3x because each corner is counted.
	;
	; When counting tetrahedra, this is 110595168 which is a 16x over-count
	; 4x because looked at each corner, and 4x because all four corners were
	; counted. So 110595168 / 16 = 6912198 (XXX or is it a 32x overcount?)
	; or is it a 32x 9x overcount?
	(define total-gene-observation-count
		(fold (lambda (gene cnt) (+ cnt (cog-count gene))) 0
			(cog-get-atoms 'GeneNode)))

	(format #t "These genes were counted a total of ~A times\n"
		total-gene-observation-count)

	; ------------------------------------------------
	; The total number of edges
	(define total-edges
		(cog-incoming-size (Predicate "gene-pair")))
	(define total-participating-edges
		(fold (lambda (pare cnt) (if (< 0 (cog-count pare)) (+ cnt 1) cnt)) 0
			(cog-incoming-set (Predicate "gene-pair"))))
	(format #t "Num participating gene-pairs ~A (out of ~A)\n"
		total-participating-edges total-edges)

	; Total count on all of the edges participating in a graph.
	; For triangles, this was 16175529, which is 9x the number of unique
	; triangles. Why 9x? A factor of 3x because each vertex is visted.
	; That is, each triangle is counted 3 times. A second factor of 3x
	; because each edge is counted once. (!) We didn't really need to count
	; each edge, by symmetry they would have shown up, eventually.
	;
	; For tetrahedra, this is 221190336.0 which is a 16x overcount (?)
	; and should be 13824396 or is this a 64x 9x overcount ???
	(define total-edge-observation-count
		(fold (lambda (pare cnt) (+ cnt (cog-count pare))) 0
			(cog-incoming-set (Predicate "gene-pair"))))

	(format #t "Gene-pairs were observed a total of ~A times\n"
		total-edge-observation-count)

	; ------------------------------------------------
	; How many triangles? We expect 1797281 of them.
	(define total-triangles
		(cog-incoming-size (Predicate "gene-triangle")))

	; But only 1701579 participate in tetrahedra
	(define total-participating-triangles
		(fold (lambda (tri cnt) (if (< 0 (cog-count tri)) (+ cnt 1) cnt)) 0
			(cog-incoming-set (Predicate "gene-triangle"))))

	(format #t "Num participating triangles ~A (out of ~A)\n"
		total-participating-triangles total-triangles)

	(define total-triangle-observation-count
		(fold (lambda (tri cnt) (+ cnt (cog-count tri))) 0
			(cog-incoming-set (Predicate "gene-triangle"))))

	(format #t "Gene-triangles were observed a total of ~A times\n"
		total-triangle-observation-count)

	; ------------------------------------------------
	; How many tetrahedra? Expect
	(define total-tetrahedra
		(cog-incoming-size (Predicate "gene-tetrahedron")))
	(format #t "Ther were ~A gene tetrahedra\n"
		total-tetrahedra)

	(define total-tetrahedra-observation-count
		(fold (lambda (pare cnt) (+ cnt (cog-count pare))) 0
			(cog-incoming-set (Predicate "gene-tetrahedron"))))

	(format #t "Gene-tetrahedra were observed a total of ~A times\n"
		total-tetrahedra-observation-count)
)

; =================================================================

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

; Genes that had positive counts.
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
; For tetrahedrons, this is 243864 (out of 365745)
;
(define participating-gene-pairs
	(filter (lambda (evlnk)
			(and (< 0 (cog-count evlnk))
				(equal? 'GeneNode (cog-type (gadr evlnk)))
				(equal? 'GeneNode (cog-type (gddr evlnk)))))
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

; When performing tirangle counting, we expect each triangle to
; be counted exactly three times, due to how the query was written.
; So .. sanity check... is this the case?
;
; How many triangles were NOT observed exactly 3 times?
; Should be zero for plain triangles.
; but something else for tetrahedra.
(length (filter (lambda (evlnk) (not (= 3 (cog-count evlnk))))
	(cog-incoming-set (Predicate "gene-triangle"))))

; ------------------------------------------------------------
; Above was for triangles. This is for tetrahedra.


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

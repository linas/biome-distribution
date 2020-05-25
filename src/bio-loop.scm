;
; bio-loop.scm
;
; Traverse genome and protemome networks, counting the number of
; various interesting gene and protein interactions. At this time,
; two graphs are examined:
;
; * Gene triangles: groups of three genes, each of which interact
;   with one-another (an order-3 clique). Takes about 110 minutes
;   on the current datasets.
; * Gene-protein pentagons: two genes that interact, and express two
;   proteins, both of which appear on the same pathway. Took 14 hours
;   on the current dataset.
;
; This started life as a modified copy of the benchmark, which in
; turn is an abstracted variant of the MOZI gene-annotation code.
;
(use-modules (srfi srfi-1))
(use-modules (opencog) (opencog exec))

; Needeed for definition of GeneNode and MoleculeNode
(use-modules (opencog bioscience))

; ------------------------------------------------------------------
; The triangle routines.

;; This defines a triangle-shaped search; one endpoint is fixed,
;; and we are looking for two other genes that interact with the
;; endpoint and form a triangle.  This is called
;; "find-output-interactors" in the MOZI code-base, and we keep
;; that name here. XXX Caution: This uses the non-symmetrized
;; edges from the original dataset. Caveat emptor!
(define (find-output-interactors gene)
	(Meet
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
		)
		(And
			(Evaluation (Predicate "interacts_with")
				(List gene (Variable "$a")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$a") (Variable "$b")))
			(Evaluation (Predicate "interacts_with")
				(List gene (Variable "$b")))
		)))

; Same as above, but using the explicitly symmetrized edges.
; It also explicitly creates the triangles as it counts.
; This helps avoid accidental double-counting and other
; hard-to-understand bugs. e.g. where the order of the distal
; gene pair is swapped.
(define (find-gene-triangles gene)
	(Query
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
		)
		(And
			(Present (Evaluation (Predicate "gene-pair")
				(Set gene (Variable "$a"))))
			(Present (Evaluation (Predicate "gene-pair")
				(Set (Variable "$a") (Variable "$b"))))
			(Present (Evaluation (Predicate "gene-pair")
				(Set (Variable "$b") gene)))
		)
		(Evaluation (Predicate "gene-triangle")
			(Set gene (Variable "$a") (Variable "$b")))
		))

;; -----------
;; Count triangles. This also creates them.
(define (count-gene-triangles gene-list)
	(define bench-secs (make-timer))
	(define batch-secs (make-timer))
	(define start-time (get-internal-real-time))
	(define ndone 0)
	(define ngen (length gene-list))
	(format #t "Begin looping over ~A genes\n" ngen)
	(for-each
		(lambda (gene)
			; Create a search pattern for each gene in the gene list.
			; (define gene (Gene gene-name))
			; (define gene-name (cog-name gene))
			; (define query (find-output-interactors gene))
			(define query (find-gene-triangles gene))

			; Perform the search
			; (define gene-secs (make-timer))
			(define result (cog-execute! query))
			(define rlist (cog-value->list result))
			(define rlen (length rlist))

			(define (mkpr a b)
				(Evaluation (Predicate "gene-pair") (Set a b)))

			; Collect up some stats. Note that this ends up
			; triple-counting everything. All counts should be
			; multiples of three! .. except for edges, which will
			; be a multiple of 9, because (3x) triple-count of triangles
			; and another (3x) cause we count every edge. So really,
			; it would have been enough to count one vertex once, and
			; one edge once. (Presuming UnorderedLink doesn't confuse
			; a single-edge counting strategy. So the below is safer.)
			(for-each
				(lambda (triangle)
					(define gene-set (gdr triangle))
					(define gene-a (cog-outgoing-atom gene-set 0))
					(define gene-b (cog-outgoing-atom gene-set 1))
					(define gene-c (cog-outgoing-atom gene-set 2))
					(define eab (mkpr gene-a gene-b))
					(define ebc (mkpr gene-b gene-c))
					(define eca (mkpr gene-c gene-a))
					(cog-inc-count! gene-a 1)
					(cog-inc-count! gene-b 1)
					(cog-inc-count! gene-c 1)
					(cog-inc-count! eab 1)
					(cog-inc-count! ebc 1)
					(cog-inc-count! eca 1)

					; Each triangle should be seen exactly 3 times.
					(cog-inc-count! triangle 1))
				rlist)

			; delete the QueryLink, too.
			(cog-delete query)
			(cog-delete-recursive (Variable "$a"))

			;; (format #t "Ran triangle ~A in ~6f seconds; got ~A results\n"
			;; 	gene-name (gene-secs) rlen)
			; (display ".")
			(set! ndone (+ ndone 1))
			(if (eq? 0 (modulo ndone 500))
				(let* ((elapsed
							(/ (- (get-internal-real-time) start-time)
								internal-time-units-per-second))
						(rate (/ ndone elapsed)))
					(format #t
						"Tri done ~A/~A in ~4f secs rate=~4f gene/sec elapsed=~6f\n"
						ndone ngen (batch-secs) rate elapsed)))
		)
		gene-list)
	(format #t "\n")
	(format #t "Finished triangle relations for ~A genes in ~8f seconds\n"
			ngen (bench-secs))

	*unspecified*
)

;; -----------
;; Explicitly create and count triangles.
;; XXX Caution: This uses the non-symmetrized edges from the original
;; dataset. Caveat emptor!
(define pointed-triangle-query
	(Query
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
			(TypedVariable (Variable "$c") (Type 'GeneNode))
		)
		(And
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$a") (Variable "$b")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$b") (Variable "$c")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$c") (Variable "$a")))
		)
		(Evaluation (Predicate "pointed_triangle")
			(List (Variable "$a") (Variable "$b") (Variable "$c")))
		))

(define (make-pointed-triangles)
	(define elapsed-secs (make-timer))
	(define pset (cog-execute! pointed-triangle-query))
	(define points (cog-value->list pset))
	(define npoints (length points))
	(format #t "Obtained ~A pointed triangles in ~6f seconds\n"
		npoints (elapsed-secs))
)

; Same as above, but not pointed; uses a set.
(define gene-triangle-query
	(Query
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
			(TypedVariable (Variable "$c") (Type 'GeneNode))
		)
		(And
			(Present (Evaluation (Predicate "gene-pair")
				(Set (Variable "$a") (Variable "$b"))))
			(Present (Evaluation (Predicate "gene-pair")
				(Set (Variable "$b") (Variable "$c"))))
			(Present (Evaluation (Predicate "gene-pair")
				(Set (Variable "$c") (Variable "$a"))))
		)
		(Evaluation (Predicate "gene-triangle")
			(Set (Variable "$a") (Variable "$b") (Variable "$c")))
		))

(define (make-gene-triangles)
	(define elapsed-secs (make-timer))
	(define tset (cog-execute! gene-triangle-query))
	(define triangles (cog-value->list tset))
	(define ntris (length triangles))
	(format #t "Obtained ~A triangles in ~6f seconds\n" ntris (elapsed-secs))
)

; =================================================================
;; -----------
;; This defines a pentagon-shaped search; one endpoint, a reaction
;; pathway, is fixed, and we are looking for two proteins that
;; participate in that pathway. These two are expressed with a pair
;; of genes that interact with one-another, forming a pentagon.
(define (pathway-gene-interactors pathway)
	(Get
		(VariableList
			(TypedVariable (Variable "$g1") (Type 'GeneNode))
			(TypedVariable (Variable "$g2") (Type 'GeneNode))
			(TypedVariable (Variable "$p1") (Type 'MoleculeNode))
			(TypedVariable (Variable "$p2") (Type 'MoleculeNode)))
		(And
			(Member (Variable "$p1") pathway)
			(Member (Variable "$p2") pathway)
			(Evaluation (Predicate "expresses") (List (Variable "$g1") (Variable "$p1")))
			(Evaluation (Predicate "expresses") (List (Variable "$g2") (Variable "$p2")))
			(Evaluation (Predicate "interacts_with") (List (Variable "$g1") (Variable "$g2")))
		)))

;; -----------
;; This defines a single edge search; one endpoint is the given
;; item, the other is a pathway.
(define (find-pathways item)
	(Get
		(TypedVariable (Variable "$p") (Type 'ConceptNode))
		(Member item (Variable "$p"))
	))

;; -----------
;; Count pentagons anchored at a pathway. A list of pathways must be supplied.
(define (count-pentagons pathways)
	(define bench-secs (make-timer))
	(define batch-secs (make-timer))
	(define start-time (get-internal-real-time))
	(define ndone 0)
	(define npath (length pathways))
	(define pr-exp (Predicate "expresses"))
	(define pr-act (Predicate "interacts_with"))

	(for-each
		(lambda (pathway)
			; Create a search pattern for each pathway in the list.
			(define query (pathway-gene-interactors pathway))

			; Perform the search
			; (define path-secs (make-timer))
			(define result (cog-execute! query))
			(define rlen (cog-arity result))

			; Collect up some stats
			(cog-inc-count! pathway rlen)
			(for-each
				(lambda (g-g-m-m)
					; The corners
					(define gene-a (cog-outgoing-atom g-g-m-m 0))
					(define gene-b (cog-outgoing-atom g-g-m-m 1))
					(define prot-a (cog-outgoing-atom g-g-m-m 2))
					(define prot-b (cog-outgoing-atom g-g-m-m 3))
					; The edges
					(define path-a (MemberLink prot-a pathway))
					(define path-b (MemberLink prot-b pathway))
					(define exp-a (Evaluation pr-exp (List gene-a prot-a)))
					(define exp-b (Evaluation pr-exp (List gene-b prot-b)))
					(define ract (Evaluation pr-act (List gene-a gene-b)))
					; Increment the counts on the corners
					(cog-inc-count! gene-a 1)
					(cog-inc-count! gene-b 1)
					(cog-inc-count! prot-a 1)
					(cog-inc-count! prot-b 1)
					; Increment the counts on the edges, too.
					(cog-inc-count! path-a 1)
					(cog-inc-count! path-b 1)
					(cog-inc-count! exp-a 1)
					(cog-inc-count! exp-b 1)
					(cog-inc-count! ract 1)
				)
				(cog-outgoing-set result))

			; delete the SetLink
			(cog-delete result)
			(cog-delete query)

			; (format #t "Ran path ~A in ~6f seconds; got ~A results\n"
			; 	(cog-name pathway) (path-secs) rlen)
			; (display ".")
			(set! ndone (+ ndone 1))
			(if (eq? 0 (modulo ndone 200))
				(let* ((elapsed
							(/ (- (get-internal-real-time) start-time)
								internal-time-units-per-second))
						(rate (/ ndone elapsed)))
					(format #t
						"Path done ~A/~A in ~4f secs rate=~4f path/sec elapsed=~6f\n"
						ndone npath (batch-secs) rate elapsed)))
		)
		pathways)
	(format #t "\n")
	(format #t "Protein expression for ~A pathways in ~6f seconds\n"
			npath (bench-secs))

	*unspecified*
)

;; -----------
(define (pathways-of-mols mol-list)
"
	Create a list of the pathways that the genes/proteins are in.
"
	(delete-dup-atoms
		(append-map
			(lambda (mol)
				(define query (find-pathways mol))
				; Perform the search
				(define path-set (cog-execute! query))
				(define pathways (cog-outgoing-set path-set))
				(cog-delete query)
				(cog-delete path-set)
				pathways
			)
			mol-list))
)

; =================================================================
; Actually do stuff. Just cut and paste from here to command line.

; (format #t "AtomSpace contents: ~A\n" (cog-report-counts))
;
; Run the triangle counting code.
; (count-triangles (cog-get-atoms 'GeneNode))
;
; (cog-incoming-size (Predicate "triangle"))
;
; Run the pentagon counting code.
; (define pathways (pathways-of-mols (cog-get-atoms 'GeneNode)))
; (length pathways)  ; 50501
; (count-pentagons pathways)  ; Takes about 16 hours...
;
; Same as above, but this time with all pathways, not just some of them.
; (define pathways (pathways-of-mols (cog-get-atoms 'MoleculeNode)))
; (length pathways)  ; 50566  so .. hardly any different ...

; ------------------------------------------------------------------

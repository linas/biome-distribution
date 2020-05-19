;
; bio-tetra.scm
;
; Traverse genome networks, counting the number of tetrahedral gene
; interactions. This is a more complex variant of `bio-loop.scm`, which
; looked at triangles and pentagons. This one explicitly looks for
; tetrahedrons.
;
(use-modules (srfi srfi-1))
(use-modules (opencog) (opencog exec))

; Needeed for definition of GeneNode and MoleculeNode
(use-modules (opencog bioscience))

; ------------------------------------------------------------------
; The tetrahedron routines.

;; This defines a tetrahedron-shaped search; one endpoint is fixed,
;; and we are looking for three other genes that interact with the
;; endpoint and form a tetrahedrone.
(define (find-gene-tetrahedron gene)
	(Get
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
			(TypedVariable (Variable "$c") (Type 'GeneNode))
		)
		(And
			(Evaluation (Predicate "interacts_with")
				(List gene (Variable "$a")))
			(Evaluation (Predicate "interacts_with")
				(List gene (Variable "$b")))
			(Evaluation (Predicate "interacts_with")
				(List gene (Variable "$c")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$a") (Variable "$b")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$b") (Variable "$c")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$c") (Variable "$a")))
		)))

;; -----------
;; Count triangles.
(define (count-triangles gene-list)
	(define bench-secs (make-timer))
	(define batch-secs (make-timer))
	(define start-time (get-internal-real-time))
	(define ndone 0)
	(define ngen (length gene-list))
	(for-each
		(lambda (gene)
			; Create a search pattern for each gene in the gene list.
			; (define gene (Gene gene-name))
			; (define gene-name (cog-name gene))
			(define query (find-gene-tetrahedron gene))

			; Perform the search
			; (define gene-secs (make-timer))
			(define result (cog-execute! query))
			(define rlen (cog-arity result))

			; Collect up some stats
			(cog-inc-count! gene rlen)
			(for-each
				(lambda (gene-triple)
					(define gene-a (cog-outgoing-atom gene-triple 0))
					(define gene-b (cog-outgoing-atom gene-triple 1))
					(define gene-c (cog-outgoing-atom gene-triple 2))
					(define act-pair (Evaluation (Predicate "interacts_with") gene-pair))

xxxxxxxxxxxxxxxxxxxx
Now what ???
					(cog-inc-count! gene-a 1)
					(cog-inc-count! gene-b 1)
					(cog-inc-count! act-pair 1))
				(cog-outgoing-set result))

			; delete the SetLink
			(cog-delete result)
			; delete the GetLink, too.
			(cog-delete query)

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
; Explicitly create and count triangles.
(define pointed-triangle-query
	(Bind
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

; Same as above, but not pointed; uses a set.
(define triangle-query
	(Bind
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
		(Evaluation (Predicate "triangle")
			(Set (Variable "$a") (Variable "$b") (Variable "$c")))
		))

(define (make-pointed-triangles)
	(define elapsed-secs (make-timer))
	(define pset (cog-execute! pointed-triangle-query))
	(define points (cog-outgoing-set pset))
	(define npoints (length points))
	(cog-delete pset)
	(format #t "Obtained ~A pointed triangles in ~6f seconds\n"
		npoints (elapsed-secs))
)

(define (make-triangles)
	(define elapsed-secs (make-timer))
	(define tset (cog-execute! triangle-query))
	(define triangles (cog-outgoing-set tset))
	(define ntris (length triangles))
	(cog-delete tset)
	(format #t "Obtained ~A triangles in ~6f seconds\n" ntris (elapsed-secs))
)

; =================================================================
; Actually do stuff. Just cut and paste from here to command line.

; (format #t "AtomSpace contents: ~A\n" (cog-report-counts))
;
; Run the triangle counting code.
; (count-triangles (cog-get-atoms 'GeneNode))

; ------------------------------------------------------------------

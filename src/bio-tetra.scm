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
;; Count tetrahedra.
(define (count-tetrahedra gene-list)
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
					(define pab (Evaluation (Predicate "interacts_with")
						(List  gene-a gene-b)))
					(define pbc (Evaluation (Predicate "interacts_with")
						(List  gene-b gene-c)))
					(define pca (Evaluation (Predicate "interacts_with")
						(List  gene-c gene-a)))

					(cog-inc-count! gene-a 1)
					(cog-inc-count! gene-b 1)
					(cog-inc-count! gene-c 1)
					(cog-inc-count! pab 1)
					(cog-inc-count! pbc 1)
					(cog-inc-count! pca 1))
				(cog-outgoing-set result))

			; delete the SetLink
			(cog-delete result)
			; delete the GetLink, too.
			(cog-delete query)

			;; (format #t "Ran tetrahedra ~A in ~6f seconds; got ~A results\n"
			;; 	gene-name (gene-secs) rlen)
			; (display ".")
			(set! ndone (+ ndone 1))
			(if (eq? 0 (modulo ndone 500))
				(let* ((elapsed
							(/ (- (get-internal-real-time) start-time)
								internal-time-units-per-second))
						(rate (/ ndone elapsed)))
					(format #t
						"Tetra done ~A/~A in ~4f secs rate=~4f gene/sec elapsed=~6f\n"
						ndone ngen (batch-secs) rate elapsed)))
		)
		gene-list)
	(format #t "\n")
	(format #t "Finished tetrahedron relations for ~A genes in ~8f seconds\n"
			ngen (bench-secs))

	*unspecified*
)

;; -----------
; Explicitly create and count tetrahedra.
(define pointed-tetrahedron-query
	(Bind
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
			(TypedVariable (Variable "$c") (Type 'GeneNode))
			(TypedVariable (Variable "$d") (Type 'GeneNode))
		)
		(And
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$a") (Variable "$b")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$b") (Variable "$c")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$c") (Variable "$a")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$a") (Variable "$d")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$b") (Variable "$d")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$c") (Variable "$d")))
		)
		(Evaluation (Predicate "pointed_tetrahedron")
			(List (Variable "$a") (Variable "$b") (Variable "$c") (Variable "$d")))
		))

; Same as above, but not pointed; uses a set.
(define tetrahedron-query
	(Bind
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
			(TypedVariable (Variable "$c") (Type 'GeneNode))
			(TypedVariable (Variable "$d") (Type 'GeneNode))
		)
		(And
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$a") (Variable "$b")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$b") (Variable "$c")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$c") (Variable "$a")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$a") (Variable "$d")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$b") (Variable "$d")))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "$c") (Variable "$d")))
		)
		(Evaluation (Predicate "tetrahedron")
			(Set (Variable "$a") (Variable "$b") (Variable "$c") (Variable "$d")))
		))

(define (make-pointed-tetrahedra)
	(define elapsed-secs (make-timer))
	(define pset (cog-execute! pointed-tetrahedron-query))
	(define points (cog-outgoing-set pset))
	(define npoints (length points))
	(cog-delete pset)
	(format #t "Obtained ~A pointed tetrahedra in ~6f seconds\n"
		npoints (elapsed-secs))
)

(define (make-tetrahedra)
	(define elapsed-secs (make-timer))
	(define tset (cog-execute! tetrahedron-query))
	(define tetrahedra (cog-outgoing-set tset))
	(define ntet (length tetrahedra))
	(cog-delete tset)
	(format #t "Obtained ~A tetrahedra in ~6f seconds\n" ntet (elapsed-secs))
)

; =================================================================
; Actually do stuff. Just cut and paste from here to command line.

; (format #t "AtomSpace contents: ~A\n" (cog-report-counts))
;
; Run the tetrahedra counting code.
; (count-tetrahedra (cog-get-atoms 'GeneNode))

; ------------------------------------------------------------------

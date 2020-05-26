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
; (use-modules (ice-9 threads))

; Needeed for definition of GeneNode and MoleculeNode
(use-modules (opencog bioscience))

; ------------------------------------------------------------------
; The tetrahedron routines.

;; This defines a tetrahedron-shaped search; one endpoint is fixed,
;; and we are looking for three other genes that interact with the
;; endpoint and form a tetrahedron.  XXX This search is not practical;
;; with the current pattern engine, this will take approx 25 cpu-days.
;; See below for an alternative.
(define (naive-find-gene-tetrahedron gene)
	(Meet
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
			(TypedVariable (Variable "$c") (Type 'GeneNode))
		)
		(And
			(Evaluation (Predicate "gene-pair")
				(Set gene (Variable "$a")))
			(Evaluation (Predicate "gene-pair")
				(Set gene (Variable "$b")))
			(Evaluation (Predicate "gene-pair")
				(Set gene (Variable "$c")))
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "$a") (Variable "$b")))
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "$b") (Variable "$c")))
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "$c") (Variable "$a")))
		)))

;; -----------
;; This defines a tetrahedron-shaped search; one endpoint is fixed,
;; and we are looking for three other genes that interact with the
;; endpoint and form a tetrahedron.  Unlike the bove search, this
;; assumes that triangles have been pre-computed.
(define (find-gene-tetrahedron gene)
	(Query
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
			(TypedVariable (Variable "$c") (Type 'GeneNode))
		)
		(Present
			(Evaluation (Predicate "gene-pair")
				(Set gene (Variable "$a")))
			(Evaluation (Predicate "gene-pair")
				(Set gene (Variable "$b")))
			(Evaluation (Predicate "gene-pair")
				(Set gene (Variable "$c")))
			(Evaluation (Predicate "gene-triangle")
				(Set (Variable "$a") (Variable "$b") (Variable "$c")))
		)
		(Evaluation (Predicate "gene-tetrahedron")
			(Set (Variable "$a") (Variable "$b") (Variable "$c") (Variable "$d")))
	))

;; -----------
;; Count tetrahedra.
(define (count-tetrahedra gene-list)
	(define bench-secs (make-timer))
	(define batch-secs (make-timer))
	(define start-time (get-internal-real-time))
	(define ndone 0)
	(define ngen (length gene-list))
	(format #t "Begin looping over ~A genes\n" ngen)
	(for-each ;; or try par-for-each
		(lambda (gene)
			; Create a search pattern for each gene in the gene list.
			; (define gene (Gene gene-name))
			; (define gene-name (cog-name gene))
			(define query (find-gene-tetrahedron gene))

			; Perform the search
			; (define gene-secs (make-timer))
			(define result (cog-execute! query))
			(define rlist (cog-value->list result))
			(define rlen (length rlist))

			(define (mkpr a b)
				(Evaluation (Predicate "gene-pair") (Set a b)))
			(define (mktri a b c)
				(Evaluation (Predicate "gene-triangle") (Set a b c)))

			; Collect up some stats. This over-counts due to
			; degeneracy, and thus wastes some CPU time. Allow
			; this as it helps provide double-checks for counting.
			(for-each
				(lambda (gene-tetra)
					(define gene-set (gdr gene-tetra))

					(define gene-a (cog-outgoing-atom gene-set 0))
					(define gene-b (cog-outgoing-atom gene-set 1))
					(define gene-c (cog-outgoing-atom gene-set 2))
					(define gene-d (cog-outgoing-atom gene-set 3))
					(define pab (mkpr gene-a gene-b))
					(define pbc (mkpr gene-b gene-c))
					(define pca (mkpr gene-c gene-a))
					(define pad (mkpr gene-a gene-d))
					(define pbd (mkpr gene-b gene-d))
					(define pcd (mkpr gene-c gene-d))

					(define tabc (mktri gene-a gene-b gene-c))
					(define tabd (mktri gene-a gene-b gene-d))
					(define tacd (mktri gene-a gene-c gene-d))
					(define tbcd (mktri gene-b gene-c gene-d))

					(cog-inc-count! gene-a 1)
					(cog-inc-count! gene-b 1)
					(cog-inc-count! gene-c 1)
					(cog-inc-count! gene-d 1)

					(cog-inc-count! pab 1)
					(cog-inc-count! pbc 1)
					(cog-inc-count! pca 1)

					(cog-inc-count! pad 1)
					(cog-inc-count! pbd 1)
					(cog-inc-count! pcd 1)

					(cog-inc-count! tabc 1)
					(cog-inc-count! tabd 1)
					(cog-inc-count! tacd 1)
					(cog-inc-count! tbcd 1)

					(cog-inc-count! gene-tetra 1))
				rlist)

			; delete the GetLink.
			(cog-delete query)
			(cog-delete-recursive (Variable "$a"))

			(set! ndone (+ ndone 1))
			(if (eq? 0 (modulo ndone 50))
				(let* ((elapsed-secs
							(/ (- (get-internal-real-time) start-time)
								internal-time-units-per-second))
						(elapsed-mins (/ elapsed-secs 60))
						(rate (/ ndone elapsed-mins)))
					(format #t
						"Tetra done ~A/~A in ~6f secs rate=~4f gene/min elapsed=~8f\n"
						ndone ngen (batch-secs) rate elapsed-secs)))
		)
		gene-list)
	(format #t "\n")
	(format #t "Finished tetrahedron relations for ~A genes in ~8f seconds\n"
			ngen (bench-secs))

	*unspecified*
)

;; -----------

;; Same as above, but without the counting.
(define (make-tetrahedra)
	(define bench-secs (make-timer))
	(define batch-secs (make-timer))
	(define start-time (get-internal-real-time))
	(define ndone 0)
	(define gene-list (cog-get-atoms 'Gene))
	(define ngen (length gene-list))
	(format #t "Begin looping over ~A genes\n" ngen)
	(for-each ;; or try par-for-each
		(lambda (gene)
			; Create a search pattern for each gene in the gene list.
			(define query (find-gene-tetrahedron gene))

			; Perform the search
			; (define gene-secs (make-timer))
			(define result (cog-execute! query))
			(define rlist (cog-value->list result))
			(define rlen (length rlist))

			; delete the GetLink.
			(cog-delete query)
			(cog-delete-recursive (Variable "$a"))

			(set! ndone (+ ndone 1))
			(if (eq? 0 (modulo ndone 50))
				(let* ((elapsed-secs
							(/ (- (get-internal-real-time) start-time)
								internal-time-units-per-second))
						(elapsed-mins (/ elapsed-secs 60))
						(rate (/ ndone elapsed-mins)))
					(format #t
						"Created Tetra ~A/~A in ~6f secs rate=~4f gene/min elapsed=~8f\n"
						ndone ngen (batch-secs) rate elapsed-secs)))
		)
		gene-list)
	(format #t "\n")
	(format #t "Finished creating tetrahedra for ~A genes in ~8f seconds\n"
			ngen (bench-secs))

	*unspecified*
)

; =================================================================
; Actually do stuff. Just cut and paste from here to command line.

; (format #t "AtomSpace contents: ~A\n" (cog-report-counts))
;
; Run the tetrahedra counting code.
; (count-tetrahedra (cog-get-atoms 'GeneNode))

; ------------------------------------------------------------------

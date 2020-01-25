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

; Performance stats timer
(define (make-timer)
	(let ((start-time (get-internal-real-time)))
		(lambda ()
			(define now (get-internal-real-time))
			(define diff (/ (- now start-time) internal-time-units-per-second))
			(set! start-time now)
			diff)))

; ------------------------------------------------------------------
; The triangle routines.

;; This defines a triangle-shaped search; one endpoint is fixed,
;; and we are looking for two other genes that interact with the
;; endpoint and form a triangle.  This is called
;; "find-output-interactors" in the MOZI code-base, and we keep
;; that name here.
(define (find-output-interactors gene)
	(Get
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
			(define query (find-output-interactors gene))

			; Perform the search
			; (define gene-secs (make-timer))
			(define result (cog-execute! query))
			(define rlen (cog-arity result))

			; Collect up some stats
			(cog-inc-count! gene rlen)
			(for-each
				(lambda (gene-pair)
					(define gene-a (cog-outgoing-atom gene-pair 0))
					(define gene-b (cog-outgoing-atom gene-pair 1))
					(define act-pair (Evaluation (Predicate "interacts_with") gene-pair))
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
		(Evaluation (Predicate "pointed_triangle")
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
;; gene, the other is a pathway.
(define (find-pathways gene)
	(Get
		(TypedVariable (Variable "$p") (Type 'ConceptNode))
		(Member gene (Variable "$p"))
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
(define (pathways-of-genes gene-list)
"
	Create a list of the pathways that the genes are in.
"
	(delete-dup-atoms
		(append-map
			(lambda (gene)
				(define query (find-pathways gene))
				; Perform the search
				(define path-set (cog-execute! query))
				(define pathways (cog-outgoing-set path-set))
				(cog-delete path-set)
				pathways
			)
			gene-list))
)

; =================================================================
; Actually do stuff. Just cut and paste from here to command line.

; (format #t "AtomSpace contents: ~A\n" (cog-report-counts))
;
; Run the triangle counting code.
; (count-triangles (cog-get-atoms 'GeneNode))

; Run the pentqagon counting code.
; (define pathways (pathways-of-genes (cog-get-atoms 'GeneNode)))
; (count-pentagons pathways)

; ------------------------------------------------------------------

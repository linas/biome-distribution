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
;; endpoint and form a triangle.
;;
;; The resulting triangles are explicitly created.
;; (This helps avoid accidental double-counting and other
;; hard-to-understand bugs. e.g. where the order of the distal
;; gene pair is swapped.)
;;
;; The original MOZI code called this `find-output-interactors`,
;; this is a fixed-up variant therof.
(define (find-gene-triangles gene)
	(Query
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
		)
		(Present
			(Evaluation (Predicate "gene-pair")
				(Set gene (Variable "$a")))
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "$a") (Variable "$b")))
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "$b") gene))
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

; =================================================================
;; -----------
;; Create triangles only, do not do any counting.
(define gene-triangle-query
	(Query
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
			(TypedVariable (Variable "$b") (Type 'GeneNode))
			(TypedVariable (Variable "$c") (Type 'GeneNode))
		)
		(Present
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "$a") (Variable "$b")))
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "$b") (Variable "$c")))
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "$c") (Variable "$a")))
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
; Actually do stuff. Just cut and paste from here to command line.

; (format #t "AtomSpace contents: ~A\n" (cog-report-counts))
;
; Run the triangle counting code.
; (count-triangles (cog-get-atoms 'GeneNode))
;
; (cog-incoming-size (Predicate "triangle"))
;
; ------------------------------------------------------------------

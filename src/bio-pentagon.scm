;
; bio-pentagon.scm
;
; Traverse genome and protemome networks, counting the number of
; various interesting gene and protein interactions. At this time,
; one graph is examined:
;
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
; Run the pentagon counting code.
; (define pathways (pathways-of-mols (cog-get-atoms 'GeneNode)))
; (length pathways)  ; 50501
; (count-pentagons pathways)  ; Takes about 16 hours...
;
; Same as above, but this time with all pathways, not just some of them.
; (define pathways (pathways-of-mols (cog-get-atoms 'MoleculeNode)))
; (length pathways)  ; 50566  so .. hardly any different ...

; ------------------------------------------------------------------

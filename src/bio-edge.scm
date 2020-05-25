;
; bio-edge.scm
;
; Traverse the genome networks, counting the number of genes that
; interact with one-another. This is the simplest part of a sequence
; of gene intraction explorations; `bio-loop.scm` looks at triangles
; (and also pentagons involving proteins) while `bio-tetra.scm` looks
; at gene-interaction tetrahedra.
;
(use-modules (srfi srfi-1))
(use-modules (opencog) (opencog exec))

; Needeed for definition of GeneNode and MoleculeNode
(use-modules (opencog bioscience))

; ------------------------------------------------------------------
; The edge-counting routines.

; This is query is designed to look like the one-dimensional version
; of the trinagle and tetrahedron query, so that all of the code can
; work the same way. Its not the most efficient way to count edges,
; but it will do.
(define (find-gene-edges gene)
	(Query
		(VariableList
			(TypedVariable (Variable "$a") (Type 'GeneNode))
		)
		(Present (Evaluation (Predicate "gene-pair")
			(Set gene (Variable "$a"))))
		(Evaluation (Predicate "gene-pair")
			(Set gene (Variable "$a")))
	))

;; -----------
;; Count the edges.
(define (count-gene-pairs gene-list)
	(define bench-secs (make-timer))
	(define batch-secs (make-timer))
	(define start-time (get-internal-real-time))
	(define ndone 0)
	(define ngen (length gene-list))
	(for-each
		(lambda (gene)
			; Create a search pattern for each gene in the gene list.
			(define query (find-gene-edges gene))

			; Perform the search
			; (define gene-secs (make-timer))
			(define result (cog-execute! query))
			(define rlist (cog-value->list result))
			(define rlen (length rlist))

			(format #t "Begin looping over %A genes\n" rlen)

			; Collect up some stats. Note that this ends up
			; double-counting everything. All counts should be
			; multiples of two!
			(for-each
				(lambda (edge)
					(define gene-set (gdr edge))
					(define gene-a (cog-outgoing-atom gene-set 0))
					(define gene-b (cog-outgoing-atom gene-set 1))
					(cog-inc-count! gene-a 1)
					(cog-inc-count! gene-b 1)

					; Each edge should be seen exactly twice.
					(cog-inc-count! edge 1))
				rlist)

			; delete the QueryLink, too.
			(cog-delete query)

			(set! ndone (+ ndone 1))
			(if (eq? 0 (modulo ndone 500))
				(let* ((elapsed
							(/ (- (get-internal-real-time) start-time)
								internal-time-units-per-second))
						(rate (/ ndone elapsed)))
					(format #t
						"Edges done ~A/~A in ~4f secs rate=~4f gene/sec elapsed=~6f\n"
						ndone ngen (batch-secs) rate elapsed)))
		)
		gene-list)
	(format #t "\n")
	(format #t "Finished counting relations for ~A genes in ~8f seconds\n"
			ngen (bench-secs))

	*unspecified*
)

; =================================================================
; Actually do stuff. Just cut and paste from here to command line.

; (format #t "AtomSpace contents: ~A\n" (cog-report-counts))
;
; Run the pair counting code.
; (count-gene-pairs (cog-get-atoms 'GeneNode))
;
; (cog-incoming-size (Predicate "gene-pair"))
;
; ------------------------------------------------------------------

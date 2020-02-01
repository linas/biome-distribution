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
   It is sorted, first.
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

; Genes that appeared in a triangular loop.
(define loop-participants
	(map (lambda (gene) (cons (cog-name gene) (cog-count gene)))
		(filter (lambda (gene) (< 0 (cog-count gene)))
			(cog-get-atoms 'GeneNode))))

(dump-to-csv loop-participants "gene-loops.csv")

; Gene pairs that appeared as edges in a triangular loop
; Some of these pairs have AnyNodes from the matrix code,
; so filter those out... 
; Some pairs have zero-count. These are interacting genes in 
; the dataset, but simply do not form a triangle.
(define gene-pairs
	(filter (lambda (evlnk)
			(and (< 0 (cog-count evlnk))
				(equal? 'GeneNode (cog-type (gadr evlnk)))
				(equal? 'GeneNode (cog-type (gddr evlnk)))))
		(cog-incoming-set (Predicate "interacts_with"))))

; Count-pairs for the gene-pairs
(define gene-pair-cnts
	(map (lambda (evelnk) (cons
		(string-concatenate
			(list (cog-name (gadr evelnk)) "-" (cog-name (gddr evelnk))))
		(cog-count evelnk)))
		gene-pairs))

(dump-to-csv gene-pair-cnts "tri-edges.csv")

; ------------------------------------------------------------
; Examination of the pentagons.

; Genes that appeared in the pentagonal loop.
(define path-genes
	(map (lambda (gene) (cons (cog-name gene) (cog-count gene)))
		(filter (lambda (gene) (< 0 (cog-count gene)))
			(cog-get-atoms 'GeneNode))))

(dump-to-csv path-genes "path-genes.csv")

(define path-proteins
	(map (lambda (protein) (cons (cog-name protein) (cog-count protein)))
		(filter (lambda (protein) (< 0 (cog-count protein)))
			(cog-get-atoms 'MoleculeNode))))

(dump-to-csv path-proteins "path-proteins.csv")

(define path-loops
	(map (lambda (pathway) (cons (cog-name pathway) (cog-count pathway)))
		(filter (lambda (pathway) (< 0 (cog-count pathway)))
			(cog-get-atoms 'ConceptNode))))

(dump-to-csv path-loops "path-loops.csv")

!# ; ---------------------------------------------------------------

; ------------------------------------------------------------------

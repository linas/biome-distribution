
;
; Notebook of cut-n-paste scripts used to generate
; the graph(s)
;
; ------------------------------------------------------------
; Gene interactions

(define gpr (Evaluation (Predicate "interacts_with") (List (Gene "FAM20C") (Gene "RNF123"))))

(define gpa (make-gene-pair-api))
(define gps (add-pair-stars gpa))
(define gpf (add-pair-freq-api gps #:nothrow #t))
(gpf 'pair-freq gpr)

(define all-gene-pairs (gps 'get-all-elts)) ; 540778 in the non-sym case,
; and 731490 in the sym case.

; Pairs with a non-zero count
(define cut-pairs
	(filter (lambda (gpr)
		(< 0 (cog-count gpr)))
	all-gene-pairs))  ; 455572 in the non-sym case, and 617530 in sym case.

; Double check that the MI value are good.
(define good-pairs
	(filter (lambda (gpr)
		(and (< 0 (cog-count gpr)) (not (inf? (gpf 'pair-fmi gpr)))))
	all-gene-pairs))  ; 455572 as above.

(gpf 'pair-fmi gpr)

; ------------------------------------------------------------
; Gene-protein expression
(load "gene-pairs.scm")
(define gea (make-expression-pair-api))
(define ges (add-pair-stars gea))

; Ick hack
(use-modules (opencog persist) (opencog persist-sql))
(sql-create "postgres:///gene_expr")
(sql-open "postgres:///gene_expr")
(batch-all-pair-mi ges)

; ---------------------------------------
; Create a "with-degeneracy" bin-count graph.
(define bins
	(bin-count good-pairs 300
		(lambda (gpr)
			(define fmi (gpf 'pair-fmi gpr))
			(if (inf? fmi) -100 fmi))
		(lambda (gpr) 1)
		-5 25))

(define fh (open-file "tri-mi.csv" "w"))
(print-bincounts-tsv bins fh)
(close fh)

; ---------------------------------------
; Create a "without-degeneracy" bin-count graph.
(define wbins
	(bin-count good-pairs 300
		(lambda (gpr) (gpf 'pair-fmi gpr))
		(lambda (gpr) (cog-count gpr))
		-5 25))

(define fh (open-file "tri-weighted-mi.csv" "w"))
(print-bincounts-tsv wbins fh)
(close fh)

; ---------------------------------------
; Make a sorted list of pairs, look at those with some minimal
; observation count, and prnt out the head of the list.

(define mi-sorted-pairs
	(sort all-gene-pairs
		(lambda (a b)
			(> (gpf 'pair-fmi a) (gpf 'pair-fmi b)))))

(define culled-pairs
	(filter (lambda (gpr)
		(and (< 4 (cog-count gpr)) (not (equal? (gadr gpr) (gddr gpr)))))
	mi-sorted-pairs))

(for-each (lambda (gpr)
	(format #t "~10A - ~10A cnt=~A   mi=~6f\n"
		(cog-name (gadr gpr))  (cog-name (gddr gpr))
		(cog-count gpr) (gpf 'pair-fmi gpr)))
	(take culled-pairs 70))

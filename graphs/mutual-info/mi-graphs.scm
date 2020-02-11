
;
; Notebook of cut-n-paste scripts used to generate
; the assorted mutual-information graphs.
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
(batch-all-pair-mi ges #f)

; ------------------------------------------------------------
; Pathway-protein pairs
(load "gene-pairs.scm")
(define ppa (make-pathway-pair-api))
(length (ppa 'all-pairs)) ; 1082860 after deleting bad chebi
(define pps (add-pair-stars ppa))
(pps 'left-basis-size) ; 14053
(pps 'right-basis-size) ; 50566
(define ppf (add-pair-freq-api pps #:nothrow #t))

; ---------------------------------------
; Create an MI bin-count graph.

(define (mi-bin-graph LLOBJ FILENAME WEIGHT-FUN)
"
  Create an MI bin-count graph for LLOBJ.
  The WEIGHT-FUN is used to assign a weight/count to each pair.
"
	(define good-pairs
		(filter (lambda (gpr)
			(and (< 0 (cog-count gpr)) (not (inf? (LLOBJ 'pair-fmi gpr)))))
			(LLOBJ 'get-all-elts)))

	(define bins
		(bin-count good-pairs 300
			(lambda (gpr)
				(define fmi (LLOBJ 'pair-fmi gpr))
				(if (inf? fmi) -100 fmi))
			WEIGHT-FUN
			-5 25))

	(define fh (open-file FILENAME "w"))
	(print-bincounts-tsv bins fh)
	(close fh)
	*unspecified*
)

; Weighted and unweighted gene-pair graphs from the triangle.
(mi-bin-graph gpf "tri-mi.csv" (lambda (pr) 1))
(mi-bin-graph gpf "tri-weighted-mi.csv" cog-count)

(mi-bin-graph ppf "pent-path-weighted-mi.csv" cog-count)

; ---------------------------------------
; Make a sorted list of pairs, look at those with some minimal
; observation count, and print out the head of the list.

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

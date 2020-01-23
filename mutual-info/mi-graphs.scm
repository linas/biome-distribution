
;
; Notebook of cut-n-paste scripts used to generate
; the graph(s)
'
(define gpr (Evaluation (Predicate "interacts_with") (List (Gene "FAM20C") (Gene "RNF123"))))

(define gpa (make-gene-pair-api))
(define gps (add-pair-stars gpa))
(define gpf (add-pair-freq-api gps #:nothrow #t))
(gpf 'pair-freq gpr)

(define all-gene-pairs (gps 'get-all-elts)) ; 540778

; Pairs with a non-zero count
(define cut-pairs
	(filter (lambda (gpr)
		(< 0 (cog-count gpr)))
	all-gene-pairs))  ; 455572

; Double check that the MI value are good.
(define good-pairs
	(filter (lambda (gpr)
		(and (< 0 (cog-count gpr)) (not (inf? (gpf 'pair-fmi gpr)))))
	all-gene-pairs))  ; 455572

(gpf 'pair-fmi gpr)

; ---------------------------------------
; Create a "with-degeneracy" bin-count graph.
(define bins
	(bin-count good-pairs 300
		(lambda (gpr)
			(define fmi (gpf 'pair-fmi gpr))
			(if (inf? fmi) -100 fmi))
		(lambda (gpr) 1)
		-10 30))

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

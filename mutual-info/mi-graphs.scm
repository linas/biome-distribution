
(define gpr (Evaluation (Predicate "interacts_with") (List (Gene "FAM20C") (Gene "RNF123"))))

(define gpa (make-gene-pair-api))
(define gps (add-pair-stars gpa))
(define gpf (add-pair-freq-api gps #:nothrow #t))
(gpf 'pair-freq gpr)

(define all-gene-pairs (gps 'get-all-elts)) ; 540778

(define cut-pairs
	(filter (lambda (gpr)
		(< 0 (cog-count gpr)))
	all-gene-pairs))  ; 455572

(define good-pairs
	(filter (lambda (gpr)
		(and (< 0 (cog-count gpr)) (not (inf? (gpf 'pair-fmi gpr)))))
	all-gene-pairs))  ; 455572




(gpf 'pair-fmi gpr)

(define bins
	(bin-count all-gene-pairs 300
		(lambda (gpr)
			(define fmi (gpf 'pair-fmi gpr))
			(if (inf? fmi) -100 fmi))
		(lambda (gpr) 1)
		-10 30))

(define fh (open-file "tri-mi.csv" "w"))
(print-bincounts-tsv bins fh)

(define mi-sorted-pairs
	(sort all-gene-pairs
		(lambda (a b)
			(> (gpf 'pair-fmi a) (gpf 'pair-fmi b)))))


(define culled-pairs
	(filter (lambda (gpr)
		(and (< 4 (cog-count gpr)) (not (equal? (gadr gpr) (gddr gpr)))))
	mi-sorted-pairs))




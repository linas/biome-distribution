
(define gpr (Evaluation (Predicate "interacts_with") (List (Gene "FAM20C") (Gene "RNF123"))))

(define gpa (make-gene-pair-api))
(define gps (add-pair-stars gpa))
(define gpf (add-pair-freq-api gps #:nothrow #t))
(gpf 'pair-freq gpr)

(define all-gene-pairs (gps 'get-all-elts)) ; 540778

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

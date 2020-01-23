; --------------------------------

(use-modules (srfi srfi-1))

; The function in question. It needs to become a runtime
; constant, for constant arguments.
(define (bar x)
	(format #t "Called bar with ~A\n" x)
	(+ x 1))

; Memoization boilerplate
(define cache (make-hash-table))
(define (int-hash INT SZ) (modulo INT SZ))
(define (int-assoc INT ILIST)
	(find (lambda (pr) (equal? INT (car pr))) ILIST))

(define-syntax foob
	(syntax-rules ()
		((foob EXP)
			(if (symbol? (quote EXP))
				(begin (display "Its a symbol\n") (bar EXP))
				(let ((junk (format #t "A hash lookup is happening for ~A\n" EXP))
						(const-val (hashx-ref int-hash int-assoc cache EXP)))
					(display "It's a constant!\n") 
      			(if (not const-val)
						(begin
							(set! const-val (bar EXP))
           				(hashx-set! int-hash int-assoc cache EXP const-val)))
					const-val)))))

; -----------------------------------------------

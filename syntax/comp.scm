

(define incr
	(let ((ctr 0))
		(lambda ()
			(set! ctr (+ 1 ctr))
			ctr)))

(define-syntax display-compile-timestamp
  (lambda (x)
    (syntax-case x ()
      ((_)
       #`(begin
          (display "The compile timestamp was: ")
          (display #,(incr))
          (newline))))))

(define (bar x)
   (format #t "Called bar with ~A\n" x)
   (+ x 1))



; (define fooj (current-time))
(define fook (display-compile-timestamp))
(define fokk (display-compile-timestamp))
(define (foom) (display-compile-timestamp))



(define (accum-time name)
"
   Maintain performance profile manually.
   Example usage:

      (define actr (accum-time \"foo\"))
      (define bctr (accum-time \"bar\"))
      (actr #:enter? #t) ; start timing foo
      (bctr #:enter? #t) ; start timing bar
      (bctr #:enter? #f) ; stop  timing bar
      (bctr #:enter? #t) ; start timing bar
      (bctr #:enter? #f) ; stop  timing bar
      (bctr #:enter? #t) ; start timing bar
      (bctr #:enter? #f) ; stop  timing bar
      (actr #:enter? #f) ; stop  timing foo
      (actr #:report? #t) ; report foo
      (bctr #:report? #t) ; report bar
"
	(let ((fname name)
			(elapsed 0)
			(calls 0)
			(start 0))
		(lambda* (#:key (enter? #f) (report? #f))
			(if report?
				(format #t "Time: ~9f secs. calls: ~A avg: ~8,1f usec/call for ~A\n"
					(* 1.0e-9 elapsed)
					calls
					(/ (* 1.0e-3 elapsed) calls)
					fname)
				(if enter?
					(set! start (get-internal-real-time))
					(begin
						(set! elapsed (+ elapsed
							(- (get-internal-real-time) start)))
						(set! calls (+ calls 1)))))))
)

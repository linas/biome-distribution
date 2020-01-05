
(define (accum-time name)
"
   Maintain performance profile manually.
   Example usage:

      (define actr (accum-time \"foo\"))
      (define bctr (accum-time \"bar\"))
      (actr #t #f) ; start timing foo
      (bctr #t #f) ; start timing bar
      (bctr #f #f) ; stop  timing bar
      (bctr #t #f) ; start timing bar
      (bctr #f #f) ; stop  timing bar
      (bctr #t #f) ; start timing bar
      (bctr #f #f) ; stop  timing bar
      (actr #f #f) ; stop  timing foo
      (actr #f #t) ; report foo
      (bctr #f #t) ; report bar
"
	(let ((fname name)
			(elapsed 0)
			(calls 0)
			(start 0))
		(lambda (enter? report?)
			(if report?
				(format #t "Time: ~6f secs. calls: ~A avg: ~2f usec/call for ~A\n"
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

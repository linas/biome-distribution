
(use-modules (ice-9 format))

(define (accum-time name)
"
   Create a performance profile with manual instrumentation.
   By installing this stub of code in various locations, performance
   timing information can be obtained for interesting code blocks.
   That is, once cna learn how much CPU time was spent in certain
   blocks of scheme code.

   Example usage:

      ; Define two performance counters.
      (define actr (accum-time \"foo\"))
      (define bctr (accum-time \"bar\"))

      ; Accumulate timing information.
      (actr #:enter? #t) ; start timing foo
      (bctr #:enter? #t) ; start timing bar
      (bctr #:enter? #f) ; stop  timing bar
      (bctr #:enter? #t) ; start timing bar
      (bctr #:enter? #f) ; stop  timing bar
      (bctr #:enter? #t) ; start timing bar
      (bctr #:enter? #f) ; stop  timing bar
      (actr #:enter? #f) ; stop  timing foo

      ; Print a performance report.
      (actr #:report? #t) ; report foo
      (bctr #:report? #t) ; report bar
"
	(let ((fname name)
			(elapsed 0)
			(calls 0)
			(start 0))
		(lambda* (#:key (enter? #f) (report? #f))
			(if report?
				(if (< 0 calls)
					(format #t
						"Time: ~9f secs. calls: ~A avg: ~8,1f usec/call for ~A\n"
						(* 1.0e-9 elapsed)
						calls
						(/ (* 1.0e-3 elapsed) calls)
						fname)
					(format #t "Zero calls to ~A\n" fname))
				(if enter?
					(set! start (get-internal-real-time))
					(begin
						(set! elapsed (+ elapsed
							(- (get-internal-real-time) start)))
						(set! calls (+ calls 1)))))))
)

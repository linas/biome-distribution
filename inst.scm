
; Maintain counts
(define (accum-time name)
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

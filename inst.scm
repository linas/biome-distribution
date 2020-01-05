
; Maintain counts
(define accum-time
	(let ((elapsed 0)
			(calls 0)
			(start 0))
		(lambda (enter? report?)
			(if report?
				(begin
				)
				(if enter?
					(set! start (get-internal-real-time))
					(begin
						(set! elapsed (+ elapsed
							(- (get-internal-real-time) start)))
						(set! calls (+ calls 1)))))))
)

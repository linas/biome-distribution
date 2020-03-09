
; Report the average time spent in GC.
(define-public report-avg-gc-cpu-time
   (let ((last-gc (gc-stats))
         (start-time (get-internal-real-time))
         (run-time (get-internal-run-time)))
      (lambda ()
         (define now (get-internal-real-time))
         (define run (get-internal-run-time))
         (define cur (gc-stats))
         (define gc-time-taken (* 1.0e-9 (- (cdar cur) (cdar last-gc))))
         (define elapsed-time (* 1.0e-9 (- now start-time)))
         (define cpu-time (* 1.0e-9 (- run run-time)))
         (define ngc (- (assoc-ref cur 'gc-times)
            (assoc-ref last-gc 'gc-times)))
         (format #t "Elapsed: ~6f secs. Rate: ~5f gc/min %cpu-GC: ~5f%  %cpu-use: ~5f%\n"
            elapsed-time
            (/ (* ngc 60) elapsed-time)
            (* 100 (/ gc-time-taken elapsed-time))
            (* 100 (/ cpu-time elapsed-time))
         )
         (set! last-gc cur)
         (set! start-time now)
         (set! run-time run))))


(set-procedure-property! report-avg-gc-cpu-time 'documentation
"
  report-avg-gc-cpu-time - Report the average time spent in GC.

  Print statistics about how much time has been spent in garbage
  collection, and how much time spent in other computational tasks.
  Resets the stats after each call, so only the stats since the
  previous call are printed.
"
)

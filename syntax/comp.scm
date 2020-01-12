
(define-syntax display-compile-timestamp
  (lambda (x)
    (syntax-case x ()
      ((_)
       #`(begin
          (display "The compile timestamp was: ")
          (display #,(current-time))
          (newline))))))


(define fooj (current-time))
(define fook (display-compile-timestamp))
(define (foom) (display-compile-timestamp))

(define (bar x)
   (format #t "Called bar with ~A\n" x)
   (+ x 1))

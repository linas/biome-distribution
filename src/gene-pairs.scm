;
; gene-pairs.scm
;
; Define a gene-pair matrix API object.
;
; Copyright (c) 2013, 2014, 2017, 2020 Linas Vepstas
;
; ---------------------------------------------------------------------
; OVERVIEW
; --------
; The objects below define API's to access natural bio-molecule pairs,
; stored in the AtomSpace, as a rank-2 matrix, i.e. as a matrix of
; (left, right) molecule-pairs.  This provides exactly the API needed
; for ; use with the `(use-modules (opencog matrix))` statistical
; analysis subsystem.
;
; To be explicit: a sparse matrix is encoded in the AtomSpace as
;
;     EvaluationLink
;         PredicateNode "some named relationship"
;         ListLink
;             MoleculeNode "some molecule or gene"
;             MoleculeNode "some other gene or protein"
;
; Although other forms are possible, as well (e.g. `InheritanceLink`s,
; etc.) Matrix entries N(x,y) are stored as counts (numbers) on the
; EvaluationLink, ; with the `x` being the first molecule, and `y`
; being the second molecule.
;
; Given the generic API, the matrx system will compute marginals N(x,*)
; and N(*,y), a grand total N(*,*) and then probabilities:
;     p(x,y) = N(x,y)/N(*,*)
; and from this, various other statistical quantities, such as mutual
; information.
;
; ---------------------------------------------------------------------
;
(use-modules (srfi srfi-1))
(use-modules (opencog))
(use-modules (opencog matrix))

; ---------------------------------------------------------------------

(define-public (make-gene-pair-api)
"
  make-gene-pair-api -- Gene-pair access methods.

  This implements a gene-pair object, where the two genes interact.
  That is, a gene pair is represented as:

    EvaluationLink
       PredicateNode \"interacts_with\"
       ListLink
          GeneNode \"SIRT1\"
          GeneNode \"LARP7\"

  After various counts, frequencies, entropies, etc pertaining to
  this particular pair are computed, they will be hung, as values,
  on the above EvaluationLink.

  The 'get-pair method returns the above EvaluationLink, if it exists.
  The 'make-pair method will create it, if it does not exist.

  In principle, the \"interacts_with\" relation is symmetric, but
  that is not how the dataset is currrently coded, so for now this
  is hacky and assymetric. XXX FIXME some day.

  Left-side counts, frequencies, etc. such as N(*,y), P(*,y) or
  log_2 P(*,y) will be placed on the following, which is returned
  by the 'left-wildcard method:

    EvaluationLink
       PredicateNode \"interacts_with\"
       ListLink
          AnyNode \"left-gene\"
          GeneNode \"LARP7\"

  The corresponding N(x,*) P(x,*) etc are hung on the atom returned
  by the 'right-wildcard method:

    EvaluationLink
       PredicateNode \"interacts_with\"
       ListLink
          GeneNode \"SIRT1\"
          AnyNode \"right-gene\"

  Finally, the 'left-type and 'right-type methods return the type
  of the the two sides of the pair.
"
	(let ((all-pairs '()))

		; Get the observational count on ATOM.
		(define (get-count ATOM) (cog-count ATOM))

		(define any-left (AnyNode "left-gene"))
		(define any-right (AnyNode "right-gene"))
		(define gene-pair-pred (PredicateNode "interacts_with"))

		(define (get-left-type) 'GeneNode)
		(define (get-right-type) 'GeneNode)
		(define (get-pair-type) 'EvaluationLink)

		; Return the atom holding the count, if it exists, else
		; return nil.
		(define (get-pair L-ATOM R-ATOM)
			(define maybe-list (cog-link 'ListLink L-ATOM R-ATOM))
			(if (null? maybe-list) '()
				(cog-link 'EvaluationLink gene-pair-pred maybe-list)))

		; Create an atom to hold the count (if it doesn't exist already).
		(define (make-pair L-ATOM R-ATOM)
			(EvaluationLink gene-pair-pred (List L-ATOM R-ATOM)))

		; Return the left member of the pair. Given the pair-atom,
		; locate the left-side atom.
		(define (get-left-element PAIR)
			(gadr PAIR))
		(define (get-right-element PAIR)
			(gddr PAIR))

		; Return the raw observational count on PAIR. If the counter for
		; PAIR does not exist (was not observed), then return 0.
		(define (get-pair-count L-ATOM R-ATOM)
			(define pr (get-pair L-ATOM R-ATOM))
			(if (null? pr) 0 (get-count pr)))

		; Caution: this unconditionally creates the wildcard pair!
		(define (get-left-wildcard WORD)
			(make-pair any-left WORD))

		; Caution: this unconditionally creates the wildcard pair!
		(define (get-right-wildcard WORD)
			(make-pair WORD any-right))

		(define (get-wild-wild)
			(make-pair any-left any-right))

		; get-all-pairs - return a list holding all of the observed
		; gene-pairs.  Caution: this can be tens of millions long!
		(define (do-get-all-pairs)
			; The list of pairs is mostly just the incoming set of the
			; predicate node. However, this does include some junk, sooo ...
			; hey, both left and right better be words.
		   (filter!
				(lambda (pair)
					(and
						(equal? 'GeneNode (cog-type (gadr pair)))
						(equal? 'GeneNode (cog-type (gddr pair)))))
				(cog-incoming-by-type gene-pair-pred 'EvaluationLink)))

		(define (get-all-pairs)
			(if (null? all-pairs) (set! all-pairs (do-get-all-pairs)))
			all-pairs)

		; fetch-gene-pairs -- fetch all counts for link-grammar
		; ANY links from the database.
		(define (fetch-gene-pairs)
			(define start-time (current-time))
			(fetch-incoming-set gene-pair-pred)
			(format #t "Elapsed time to load gene pairs: ~A secs\n"
				(- (current-time) start-time))
		)

		; Delete the pairs from the atomspace AND the database.
		; But only those that are currently in the atomspace are
		; deleted; if any are hiding in the database, they will not be
		; touched.
		(define (delete-gene-pairs)
			(define start-time (current-time))
			(for-each (lambda (PAIR) (cog-delete-recursive (gdr PAIR)))
				(cog-incoming-set gene-pair-pred))
			(cog-delete gene-pair-pred)
			(cog-delete any-left)
			(cog-delete any-right)
			(format #t "Elapsed time to delete gene pairs: ~A secs\n"
				(- (current-time) start-time))
		)

		; Methods on the object
		(lambda (message . args)
			(apply (case message
					((name) (lambda () "Gene interacts_with Pairs"))
					((id)   (lambda () "Gene-pairs"))
					((left-type) get-left-type)
					((right-type) get-right-type)
					((pair-type) get-pair-type)
					((pair-count) get-pair-count)
					((get-pair) get-pair)
					((get-count) get-count)
					((make-pair) make-pair)
					((left-element) get-left-element)
					((right-element) get-right-element)
					((left-wildcard) get-left-wildcard)
					((right-wildcard) get-right-wildcard)
					((wild-wild) get-wild-wild)
					((all-pairs) get-all-pairs)
					((fetch-pairs) fetch-gene-pairs)
					((delete-pairs) delete-gene-pairs)
					((provides) (lambda (symb) #f))
					((filters?) (lambda () #f))
					(else (error "Bad method call on gene-pair-api:" message)))
				args)))
)


; ---------------------------------------------------------------------
; Handy-dandy main entry points. These compute mutual information.
;
; XXX FIXME - these will fail if an SQL database is not open.
; For now, the work-around is to create a database, and open it.
; As follows:
;
; (use-modules (opencog persist) (opencog persist-sql))
; (sql-create "postgres:///gene_pairs")
; (sql-open "postgres:///gene_pairs")
;

(define-public (batch-pairs LLOBJ)
	(cog-report-counts)
	(batch-all-pair-mi LLOBJ)
	(print-matrix-summary-report LLOBJ)
)

(define-public (batch-gene-pairs)
	(batch-pairs (make-gene-pair-api))
)

; ---------------------------------------------------------------------
;
; Debugging notes and cheat-sheet.
; (define gpr (Evaluation (Predicate "interacts_with") (List (Gene "FAM20C") (Gene "RNF123"))))
; (define gpa (make-gene-pair-api))
; (define gps (add-pair-stars gpa))
; (define gpf (add-pair-freq-api gps #:nothrow #t))
; (gpf 'pair-freq gpr)
;
; (define all-gene-pairs (gps 'get-all-elts)) ; 540778
; (define any-left (gps 'left-wild-pairs))    ; 18074
; (define any-right (gps 'right-wild-pairs))  ; 19648
; (define evs (cog-incoming-set (Predicate "interacts_with"))) ; 578501
; (define hmm (atoms-subtract evs all-gene-pairs))
; (define h2 (atoms-subtract hmm any-left))
; (define wtf (atoms-subtract h2 any-right))

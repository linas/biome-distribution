;
; gene-pairs.scm
;
; Define gene-pair access API objects.
; Batch-compute the mutual information of pairs of genes.
;
; Copyright (c) 2013, 2014, 2017, 2020 Linas Vepstas
;
; ---------------------------------------------------------------------
; OVERVIEW
; --------
; The objects below define API's to access natural bio-molecule pairs,
; stored in the AtomSpace, as a rank-2 matrix, i.e. as a matrix of
; (left, right) molecule-pairs.  This provides exactly the API needed for
; use with the `(use-modules (opencog matrix))` statistical analysis
; subsystem.
;
; Given a generic API, the `(opencog matrix)` can do things such as
; computing the Yuret-style lexical attraction between pairs of molecules.
; (See `compute-mi.scm` for more detail about what is computed, and how.)
;
; Given the generic API, there is a handful of small scripts, at the
; bottom of this file, that will perform the MI calculations as a batch
; job.  As a batch job, and may take hours to complete. The results are
; stored in the currently-open database, for future reference.
;
; One structure, among several, in which the pair counts are held,
; is of the form
;
;     EvaluationLink
;         LinkGrammarRelationshipNode "ANY"
;         ListLink
;             WordNode "some-word"
;             WordNode "other-word"
;
; After they've been computed, the values for N(w,*) and N(*,w) can be
; fetched with the `get-left-count-str` and `get-right-count-str`
; routines, below.  The value for N(*,*) can be gotten by calling
; `total-pair-observations`.
;
; The counting done in `link-pipeline.scm` keeps track of several
; different types of pair information.  Besides the above, it also
; counts these things:
;
;     EvaluationLink
;         PredicateNode "*-Sentence Word Pair-*"
;         ListLink
;             WordNode "lefty"
;             WordNode "righty"
;
;     ExecutionLink
;         SchemaNode "*-Pair Distance-*"
;         ListLink
;             WordNode "lefty"
;             WordNode "righty"
;         NumberNode 3
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

  Left-side counts, frequencies, etc. such as N(*,y) P(*,y) or
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
					(else (error "Bad method call on ANY-link:" message)))
				args)))
)


; ---------------------------------------------------------------------
; Handy-dandy main entry points.

(define-public (batch-pairs LLOBJ)
	(cog-report-counts)
	(batch-all-pair-mi LLOBJ)
	(print-matrix-summary-report LLOBJ)
)

(define-public (batch-gene-pairs)
	(batch-pairs (make-gene-pair-api))
)

; ---------------------------------------------------------------------

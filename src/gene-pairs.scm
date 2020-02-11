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
; for use with the `(use-modules (opencog matrix))` statistical
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
; EvaluationLink, with the `x` being the first molecule, and `y`
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

  This implements a matrix object representing gene-pairs, denoting
  that the two genes interact.  A gene pair is represented as:

    EvaluationLink
       PredicateNode \"interacts_with\"
       ListLink
          GeneNode \"SIRT1\"
          GeneNode \"LARP7\"

  This Atom (the EvaluationLink) be used to record counts, frequencies,
  entropies, etc pertaining to this particular pair, simply by placing
  them, as values, on the above EvaluationLink.

  The 'get-pair method returns the above EvaluationLink, if it exists.
  The 'make-pair method will create it, if it does not exist.

  In principle, the \"interacts_with\" relation is symmetric, but
  that is not how the dataset is currrently coded, so for now this
  is hacky and asymmetric. XXX FIXME some day.

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
"
	(make-evaluation-pair-api
		(PredicateNode "interacts_with")
		'GeneNode 'MoleculeNode
		(AnyNode "left-gene") (AnyNode "right-gene")
		"Gene-pairs"
		"Gene pairs, predicate `interacts_with`")
)

; ---------------------------------------------------------------------

(define-public (make-expression-pair-api)
"
  make-expression-pair-api -- Gene-protein expression access methods.

  This implements a matrix object, connecting genes and proteins.
  The rows of the matrix are genes, the columns are proteins.
  These are represented as:

    EvaluationLink
       PredicateNode \"expresses\"
       ListLink
          MoleculeNode \"Uniprot:E9PC49\"
          GeneNode \"SIRT1\"

  This Atom (the EvaluationLink) be used to record counts, frequencies,
  entropies, etc pertaining to this particular pair, simply by placing
  them, as values, on the above EvaluationLink.

  The 'get-pair method returns the above EvaluationLink, if it exists.
  The 'make-pair method will create it, if it does not exist.

  Left-side counts, frequencies, etc. such as N(*,y), P(*,y) or
  log_2 P(*,y) will be placed on the following, which is returned
  by the 'left-wildcard method:

    EvaluationLink
       PredicateNode \"expresses\"
       ListLink
          AnyNode \"left-protein\"
          GeneNode \"SIRT1\"

  The corresponding N(x,*) P(x,*) etc are hung on the atom returned
  by the 'right-wildcard method:

    EvaluationLink
       PredicateNode \"expresses\"
       ListLink
          MoleculeNode \"Uniprot:E9PC49\"
          AnyNode \"right-gene\"
"
	(make-evaluation-pair-api
		(PredicateNode "expresses")
		'GeneNode 'MoleculeNode
		(AnyNode "left-gene") (AnyNode "right-protein")
		"Gene-expression"
		"Gene-Protein pairs, predicate `expresses`")
)

; ---------------------------------------------------------------------

(define-public (make-pathway-pair-api)
"
  make-pathway-pair-api -- Pathway-protein pair access methods.

  This implements a matrix object, connecting pathways and proteins.
  The rows of the matrix are pathways, the columns are proteins.
  These are represented as:

    MemberLink
       MoleculeNode \"ChEBI:16908\"
       ConceptNode \"SMP0027226\"

  This Atom (the MemberLink) be used to record counts, frequencies,
  entropies, etc pertaining to this particular pair, simply by placing
  them, as values, on the MemberLink.

  The 'get-pair method returns the above MemberLink, if it exists.
  The 'make-pair method will create it, if it does not exist.

  Left-side counts, frequencies, etc. such as N(*,y), P(*,y) or
  log_2 P(*,y) will be placed on the following, which is returned
  by the 'left-wildcard method:

    MemberLink
       AnyNode \"left-protein\"
       ConceptNode \"SMP0027226\"

  The corresponding N(x,*) P(x,*) etc are hung on the atom returned
  by the 'right-wildcard method:

    MemberLink
      MoleculeNode \"ChEBI:16908\"
      AnyNode \"right-pathway\"
"
	(let ((all-pairs '()))

		; Get the observational count on ATOM.
		(define (get-count ATOM) (cog-count ATOM))

		(define any-left (AnyNode "left-protein"))
		(define any-right (AnyNode "right-pathway"))

		(define (get-left-type) 'MoleculeNode)
		(define (get-right-type) 'ConceptNode)
		(define (get-pair-type) 'MemberLink)

		; Return the atom holding the count, if it exists, else
		; return nil.
		(define (get-pair L-ATOM R-ATOM)
			(cog-link 'MemberLink L-ATOM R-ATOM))

		; Create an atom to hold the count (if it doesn't exist already).
		(define (make-pair L-ATOM R-ATOM) (MemberLink L-ATOM R-ATOM))

		; Return the left member of the pair. Given the pair-atom,
		; locate the left-side atom.
		(define (get-left-element PAIR) (gar PAIR))
		(define (get-right-element PAIR) (gdr PAIR))

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

		; get-all-pairs - return a list holding all of the
		; pathway-protein pairs.
		(define (do-get-all-pairs)
		   (filter!
				(lambda (pair)
					(define path (gdr pair))
					(and
						(equal? 'MoleculeNode (cog-type (gar pair)))
						(equal? 'ConceptNode (cog-type path))
						(or
							(string-prefix? "SMP" (cog-name path))
							(string-prefix? "R-HSA-"(cog-name path)))))
				(cog-get-atoms 'MemberLink)))

		(define (get-all-pairs)
			(if (null? all-pairs) (set! all-pairs (do-get-all-pairs)))
			all-pairs)

		; fetch-gene-pairs -- fetch all counts from the database.
		; NOT IMPLEMENTED
		(define (fetch-pairs) #f)

		; Delete the pathway-protein pairs from the atomspace AND the database.
		; But only those that are currently in the atomspace are
		; deleted; if any are hiding in the database, they will not be
		; touched.
		(define (delete-pairs)
			(define start-time (current-time))
			(for-each (lambda (PAIR) (cog-delete-recursive PAIR))
				(get-all-pairs))
			(cog-delete any-left)
			(cog-delete any-right)
			(format #t "Elapsed time to delete pathway-proteins: ~A secs\n"
				(- (current-time) start-time))
		)

		; Methods on the object
		(lambda (message . args)
			(apply (case message
					((name) (lambda () "Pathway-protein pairs, MemberLink"))
					((id)   (lambda () "pathway-protein-pairs"))
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
					((fetch-pairs) fetch-pairs)
					((delete-pairs) delete-pairs)
					((provides) (lambda (symb) #f))
					((filters?) (lambda () #f))
					(else (error "Bad method call on pathway-pair-api:" message)))
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

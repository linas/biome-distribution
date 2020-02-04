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

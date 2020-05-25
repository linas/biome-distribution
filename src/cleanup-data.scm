;
; cleanup-data.scm
;
; The raw data is a tad messy and problematic. The functions here fix
; assorted issues.
;
(use-modules (opencog exec))
(use-modules (srfi srfi-1))

; --------------------
; Performance stats timer
(define (make-timer)
	(let ((start-time (get-internal-real-time)))
		(lambda ()
			(define now (get-internal-real-time))
			(define diff (/ (- now start-time) internal-time-units-per-second))
			(set! start-time now)
			diff)))

; --------------------
(define (delete-go-nodes)
"
  Delete the GO (GeneOnotology) nodes, as they are not pathways.
  They way they are encoded interferes with pathway nodes.
"
	(for-each
		(lambda (cpt)
			(if (string-contains (cog-name cpt) "GO:")
				(cog-delete-recursive cpt)))
		(cog-get-atoms 'ConceptNode)))

; --------------------
(define (delete-self-interaction)
"
  Many genes are marked as interacting with themselves.
  Delete these, they screw up the topology of the searches.
"
	(define selfie-q
		(Get (List (Variable "$x") (Variable "$x"))))
	(define selfie-set (cog-execute! selfie-q))
	(define selfies (cog-outgoing-set selfie-set))
	(cog-delete selfie-set)
	(format #t "Deleting ~A self-interactions\n" (length selfies))
	(for-each
		(lambda (gene) (cog-delete-recursive (List gene gene)))
		selfies))

; --------------------
(define (delete-bad-chebi)
"
  Delete (MoleculeNode \"ChEBI:nan\") and all links that contain it.
  This is not a valid protein. Also delete some other junk.
"
	(cog-delete-recursive (MoleculeNode "ChEBI:nan"))
	(cog-delete-recursive (ConceptNode "SMPD1 "))
	(cog-delete-recursive (ConceptNode "SMPD2 "))
	(cog-delete-recursive (ConceptNode "SMPD3 "))
	(cog-delete-recursive (ConceptNode "SMPD4 "))

	*unspecified*
)

; --------------------
(define (count-gene-interactions)
"
  Count the number of symmetric and non-symmetric gene-pair interactions
  in the dataset.  The gene interactions use the asymmetric ListLink to
  denote gene pairs. This counts the total number of interactions, the
  number of paired and unpaired interactions, and the number of
  self-interactions.
"

	; Total number of interactions
	(define count-q
		(Get
			(VariableList
				(TypedVariable (Variable "g1") (Type 'GeneNode))
				(TypedVariable (Variable "g2") (Type 'GeneNode)))
			(Evaluation (Predicate "interacts_with")
					(List (Variable "g1") (Variable "g2")))))

	; Symmetric interactions
	(define sym-q
		(Get
			(VariableList
				(TypedVariable (Variable "g1") (Type 'GeneNode))
				(TypedVariable (Variable "g2") (Type 'GeneNode)))
			(And
				(Evaluation (Predicate "interacts_with")
						(List (Variable "g1") (Variable "g2")))
				(Evaluation (Predicate "interacts_with")
						(List (Variable "g2") (Variable "g1"))))))

	; Self-interactions
	(define self-q
		(Get
			(TypedVariable (Variable "g1") (Type 'GeneNode))
			(Evaluation (Predicate "interacts_with")
					(List (Variable "g1") (Variable "g1")))))

	(define cnt-set (cog-execute! count-q))
	(define n-acts (length (cog-outgoing-set cnt-set)))
	(cog-delete cnt-set)

	(format #t "Found ~A gene interactions\n" n-acts)

	(define sym-set (cog-execute! sym-q))
	(define n-sym (length (cog-outgoing-set sym-set)))
	(cog-delete sym-set)

	(format #t "Found ~A symmetric (paired) gene interactions\n" n-sym)

	(define self-set (cog-execute! self-q))
	(define n-self (length (cog-outgoing-set self-set)))
	(cog-delete self-set)
	(format #t "Found ~A self-interacting  genes\n" n-self)

	(cog-delete-recursive (Variable "g1"))
	(cog-delete-recursive (Variable "g2"))

	(define n-asym (- n-acts n-sym))
	(format #t "Conclude: there are ~A asymmetric interctions\n" n-asym)

	(define n-edge (/ (- (+ n-acts n-asym) n-self) 2))
	(format #t "Conclude: there are ~A symmetrized interactions\n" n-edge)
	*unspecified*
)

; --------------------
(define (symmetrize-gene-interactions)
"
  The gene interactions use the asymmetric ListLink to denote
  gene pairs. But gene interactions are (meant to be) symmetric, so
  they should have used the SetLink. Oh well. At this time, it is
  convenient to use the ListLink during pattern searches; but in this
  case, the interactions should be symmetrized. That's what this does.
"
	(define interact-q
		(Query
			(VariableList
				(TypedVariable (Variable "g1") (Type 'GeneNode))
				(TypedVariable (Variable "g2") (Type 'GeneNode)))
			(Present
				(Evaluation (Predicate "interacts_with")
					(List (Variable "g1") (Variable "g2"))))
			(Evaluation (Predicate "interacts_with")
				(List (Variable "g2") (Variable "g1")))))

	(define sym-set (cog-execute! interact-q))
	(format #t "Found ~A gene interactions\n"
		(length (cog-value->list sym-set)))
	(cog-delete interact-q)
	*unspecified*
)

; --------------------
(define (make-gene-pairs)
"
  The gene interactions use the asymmetric ListLink to denote
  gene pairs. But gene interactions are (meant to be) symmetric, so
  they should have used the SetLink. We make this now.

  See also: make-triangles, make-tetrahedra
"
	(define make-sym-pairs
		(Query
			(VariableList
				(TypedVariable (Variable "g1") (Type 'GeneNode))
				(TypedVariable (Variable "g2") (Type 'GeneNode)))
			(Present
				(Evaluation (Predicate "interacts_with")
					(List (Variable "g1") (Variable "g2"))))
			(Evaluation (Predicate "gene-pair")
				(Set (Variable "g1") (Variable "g2")))))

	(define elapsed-secs (make-timer))
	(define sym-set (cog-execute! make-sym-pairs))
	(format #t "Created ~A symmetric gene-pairs in ~6f seconds\n"
		(length (cog-value->list sym-set)) (elapsed-secs))
	(cog-delete make-sym-pairs)
	*unspecified*
)

; --------------------
(define (delete-all-but-interactions)
"
  Remove pretty much everything that is not a gene or a protein
  interacting with one-another.
  The goal is to reduce the atomspace to something manageable and more
  responsive during gene graph searches.
  This assumes that `(make-gene-pairs)` has alreay run.
"
	; What about the first two?
	(cog-delete-recursive (PredicateNode "transcribed_to"))
	(cog-delete-recursive (PredicateNode "translated_to"))
	(cog-delete-recursive (PredicateNode "has_location"))
	(cog-delete-recursive (PredicateNode "has_name"))
	(cog-delete-recursive (PredicateNode "has_pubmedID"))
	(cog-delete-recursive (PredicateNode "GO_name"))
	(cog-delete-recursive (PredicateNode "has_biogridID"))
	(cog-delete-recursive (PredicateNode "GO_namespace"))
	(cog-delete-recursive (PredicateNode "has_entrez_id"))

	; The symmetrized gene-pair has rendered this useless.
	(cog-delete-recursive (PredicateNode "interacts_with"))

	; The pentagons depend on interacting genes, so kill all genes
	; that aren't in some gene-pair.
	(define all-genes (cog-get-atoms 'GeneNode))
	(define interacting-genes
		(filter
			(lambda (gene) (< 0 (cog-incoming-size-by-type gene 'Set)))
			all-genes))

	(for-each cog-delete-recursive
		(atoms-subtract all-genes interacting-genes))

	; The above will orphan many ListLinks. Delete them.
	(for-each cog-delete
		(filter
			(lambda (lst) (= 0 (cog-incoming-size lst)))
			(cog-get-atoms 'List)))

	; We are not looking at InheritanceLinks for anything
	(for-each
		(lambda (misc)
			(if (= 0 (cog-incoming-size misc)) (cog-delete misc)))
		(cog-get-atoms 'Inheritance))

	; Once the ListLinks are gone, then orphan nodes show up.
	(for-each
		(lambda (misc)
			(if (= 0 (cog-incoming-size misc)) (cog-delete misc)))
		(cog-get-atoms 'Concept))

	(for-each
		(lambda (misc)
			(if (= 0 (cog-incoming-size misc)) (cog-delete misc)))
		(cog-get-atoms 'Molecule))

	(for-each
		(lambda (gene)
			(if (= 0 (cog-incoming-size gene)) (cog-delete gene)))
		(cog-get-atoms 'Gene))
)

; --------------------
(define (delete-all-but-gene-interactions)
"
  Remove pretty much everything that is not a gene that belongs to
  a gene-pair.  The goal is to reduce the atomspace to something
  manageable and more responsive during gene graph searches.
  This assumes that `(make-gene-pairs)` has alreay run.
"
	(define elapsed-secs (make-timer))
	(for-each cog-delete-recursive (cog-get-atoms 'Molecule))
	(for-each cog-delete-recursive (cog-get-atoms 'Concept))
	(for-each cog-delete-recursive (cog-get-atoms 'List))

	(for-each
		(lambda (pred)
			(if (not (equal? pred (Predicate "gene-pair")))
				(cog-delete-recursive pred)))
		(cog-get-atoms 'Predicate))

	(for-each
		(lambda (gene)
			(if (= 0 (cog-incoming-size gene)) (cog-delete gene)))
		(cog-get-atoms 'Gene))

	(format #t "Cleaned out non-genomic data ~6f seconds\n" (elapsed-secs))
	(format #t "What's left: ~A\n" (cog-report-counts))
	*unspecified*
)

; --------------------
(define (delete-simple-tv)
"
  Delete the SimpleTruthValues on all atoms in the atomspace.
  The problem is that calling `get-count` on a SimpleTruthValue
  returns garbage, thus messing up statistics. Unfortunately,
  this cannot be fixed, because the PLN book documents the garbage;
  its part of the spec. Whoops.
"
	(define elapsed-secs (make-timer))
	; Setting to (stv 1 0) sets it to DEFAULT_TV, which frees
	; the RAM in the AtomSpace.
	(for-each
		(lambda (ATOM)
			(if (not (cog-ctv? (cog-tv ATOM)))
				(cog-set-tv! ATOM (stv 1 0))))
		(cog-get-atoms 'Atom #t))
	(format #t "Removed SimpleTV in ~6f seconds\n" (elapsed-secs))
	*unspecified*
)

; --------------------

;
; cleanup-data.scm
;
; The raw data is a tad messy and problematic. The functions here fix
; assorted issues.
;
(use-modules (opencog exec))

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

(define (symmetrize-gene-interactions)
"
  The gene interactions use the asymmetric ListLink to denote
  gene pairs. But gene interactions are symmetrix, so they should
  have used the SetLink. Oh well. At this time, it is convenient
  to use the ListLink during pattern searches; but in this case,
  the interactions should be symmetrized. That's what this does.
"
	(define interact-q
		(Bind
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
		(length (cog-outgoing-set sym-set)))
	(cog-delete sym-set)
)

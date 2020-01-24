;
; cleanup-data.scm
;
; The raw data is a tad messy and problematic. The functions here fix
; assorted issues.
;
(use-modules (opencog exec))

; Delete the GO (GeneOnotology) nodes, as they are not pathways.
; They way they are encoded interferes with pathway nodes.
(define (delete-go-nodes)
	(for-each
		(lambda (cpt)
			(if (string-contains (cog-name cpt) "GO:")
				(cog-delete-recursive cpt)))
		(cog-get-atoms 'ConceptNode)))

; Many genes are marked as interacting with themselves.
; Delete these, the screw up the topology of the searches.
(define (delete-self-interaction)
	(define selfie-q
		(Get (List (Variable "$x") (Variable "$x"))))
	(define selfie-set (cog-execute! selfie-q))
	(define selfies (cog-outgoing-set selfie-set))
	(cog-delete selfie-set)
	(for-each
		(lambda (gene) (cog-delete-recursive (List gene gene)))
		selfies))

;
; cleanup-data.scm
;
; The raw data is a tad messy and problematic. The functions here fix
; assorted issues.
;

; Delete the GO (GeneOnotology) nodes, as they are not pathways.
; They way they are encoded interferes with pathway nodes.
(define (delete-go-nodes)
	(for-each
		(lambda (cpt)
			(if (string-contains (cog-name cpt) "GO:")
				(cog-delete-recursive cpt)))
		(cog-get-atoms 'ConceptNode)))

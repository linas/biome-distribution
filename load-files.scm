
(use-modules (opencog))
; (use-modules (annotation))
(use-modules (opencog bioscience))
(use-modules (opencog cogserver))

(define (loaf fn)
	(define start (current-time))
	(define ssz (count-all))
	(define path (string-append "/home/ubuntu/datasets/" fn))
	(define jnk (primitive-load path))
	(define delapse (- (current-time) start))
	(define elapse (if (eq? delapse 0) 1 delapse))
	(define delta (- (count-all) ssz))
	(define rate (/ (exact->inexact delta) elapse))
	(format #t "Took ~A secs to load ~A\n" elapse fn)
	(format #t "Loaded ~A atoms (~A per sec) total atoms=~A\n"
		delta rate (count-all))
	*unspecified*
)

(use-modules (ice-9 threads))

(define file-list (list
 "biogridgene2uniprot.scm"
 "biogrid_gene_gene_174.scm"
 "ChEBI2Reactome_PE_Pathway.txt.scm"
 "current_symbols.scm"
 "entrez_to_protein.scm"
 "GO_annotation.scm"
 "GO_without_def.scm"
 "NCBI2Reactome_PE_Pathway.txt.scm"
 "reactome.scm"
 "smpdb_chebi_wname.scm"
 "smpdb_protein.scm"
 "uniprot2GO.scm"
 "UniProt2Reactome_PE_Pathway.txt.scm"
))

(define current-2019-list (list
 "current-2019-12-31/biogridgene2uniprot.scm"
 "current-2019-12-31/biogrid_gene_gene_3.5.177.scm"
 "current-2019-12-31/ChEBI2Reactome_PE_Pathway.txt.scm"
 "current-2019-12-31/entrez_to_protein.scm"
 "current-2019-12-31/GO.scm"
 "current-2019-12-31/GO_annotation.scm"
 "current-2019-12-31/NCBI2Reactome_PE_Pathway.txt.scm"
 "current-2019-12-31/reactome.scm"
 "current-2019-12-31/uniprot2GO.scm"
 "current-2019-12-31/UniProt2Reactome_PE_Pathway.txt.scm"
))

(define current-list (list
	"current/ChEBI2Reactome_PE_Pathway.txt.scm"
	"current/biogrid_gene_gene_3.5.177.scm"
	"current/noncodingRNA.scm"
	"current/smpdb_protein.scm"
	"current/NCBI2Reactome_PE_Pathway.txt.scm"
	"current/uniprot2GO.scm"
	"current/UniProt2Reactome_PE_Pathway.txt.scm"
	"current/reactome.scm"
	"current/GO_annotation.scm"
	"current/entrez_to_protein.scm"
	"current/current_symbols.scm"
	"current/smpdb_chebi_wname.scm"
	"current/biogridgene2uniprot.scm"
	"current/codingRNA.scm"
	"current/GO.scm"
))

(define (serial-load-all)
	(define start (current-time))
	; Using more than 2 threads results in net de-acceleration!
	; (n-par-for-each 2 loaf file-list)
	(for-each loaf current-list)
	(format #t "\nLoaded all the files in ~A seconds\n" (- (current-time) start))
	#f
)

(define (load-all)
	(define start (current-time))
	; Using more than 2 threads results in net de-acceleration!
	; (n-par-for-each 2 loaf file-list)
	(n-par-for-each 2 loaf current-list)
	(format #t "\nLoaded all the files in ~A seconds\n" (- (current-time) start))
	#f
)

(start-cogserver)
*unspecified*

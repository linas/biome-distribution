
Mutual information
------------------

The `gene-pairs.scm` provides the rank-2 matrix API on top of the
atomspace. Currently supported are:
```
   Evaluation (count=n)
      Predicate "interacts_with"
		List
         Gene "FOO"
         Gene "BAR"
```
where `(count=n)` is a CountTrutheValue holding observation counts for
the gene-pair (FOO, BAR). This is a matrix, of course, with rows FOO
and columns BAR. The matrix API is provided by `make-gene-pair-api`.

Similarly, the `make-expression-pair-api` wrapper provides access to
gene-protein expressions recorded as
```
   Evaluation (count=n)
      Predicate "expresses"
		List
         Gene "FOO"
         Molecule "BAR"
```

The `make-pathway-pair-api` wrapper provides access to pathway-protein
pair recorded as
```
   Member (count=n)
      Molecule "Protein:FOO"
      Concept "pathway-BAR"
```

Mutual information is computed by running `batch-all-pair-mi` for each
of these matrix objects.

The file `mi-graphs.scm` is some hacky cut-n-paste scriptlets used to
generate graphs.

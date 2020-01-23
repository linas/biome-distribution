
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
and columns BAR.


The file `mi-graphs.scm` is some hacky cut-n-paste scriptlets used to
generate graphs.

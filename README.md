
Biome Distributions
-------------------
Research concerning the statistical distribution of genetic
interactions, the proteins expressed by genes, and their
participation in reaction pathways.

The primary result is a paper, being prepared for publication.
It is in the [paper](./paper) directory. The
[PDF is here](./paper/biome-distributions.pdf).

It was naively hypothesized that genome/proteome reaction pathways
form a scale-free network, and thus would have a Zipfian distribution.
Much to our surprise, this is not the case! It seems like *everything*
follows a square-root Zipfian distribution! I do not know of any
network theory or biology theory that would explain this, so it is
a surprise.

An exploration of the mutual information of interaction pathways is
also performed. It appears that these are easily fit with a bimodal
Gaussian distribution.

This is for human genome/reactome data. I don't doubt that the results
are generic in biology.

### Directory layout
Other directories here:

* [diary](./diary) A research diary of notes.
* [graphs](./graphs) Graphs and the tools used to prepare them.
* [src](./src) Source code for loading, processing and analyzing
  the genomic and proteomic data, including reactome data, etc.

### License
<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.

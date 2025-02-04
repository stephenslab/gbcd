# gbcd

R package implementing the Generalized Binary Covariance Decomposition
(GBCD) method for dissecting tumor transcriptional heterogeneity
from multi-tumor single-cell RNA-sequencing data.
	
If you find a bug, or you have a question or feedback on this software,
please post an [issue][issues].

## Citing this work

If you find the gbcd R package or any of the source code in this
repository useful for your work, please cite:

> Y. Liu, P. Carbonetto, J. Willwerscheid, S. A. Oakes, K. F. Macleod
> and M. Stephens. [Dissecting tumor transcriptional
> heterogeneity from single-cell RNA-seq data by generalized binary
> covariance decomposition.](https://doi.org/10.1038/s41588-024-01997-z)
> *Nature Genetics* **57**, 263-273 (2025).

## Quick Start

Install and load the package from GitHub: 

```r
# install.packages("remotes")
remotes::install_github("stephenslab/gbcd")
```

See the [package documentation][pkgdown] and the [vignette][vignette]
for more information.

[pkgdown]: https://stephenslab.github.io/gbcd/
[vignette]: https://stephenslab.github.io/gbcd/articles/hnscc.html
[gbcd-biorxiv]: https://doi.org/10.1101/2023.08.15.553436
[issues]: https://github.com/stephenslab/gbcd/issues

[![Build Status](https://travis-ci.com/fungenomics/cytobox.svg?branch=master)](https://travis-ci.com/fungenomics/cytobox)

# cytobox: README
A toolkit for analyzing single cell RNA-seq data created by the Kleinman Lab. The legacy package,
which is no longer being developed, is available at https://github.com/sjessa/cytobox/

## Installation

You will need to install [`devtools`](https://cran.r-project.org/web/packages/devtools/), and then run:

```r
devtools::install_github("fungenomics/cytobox", build_vignettes = TRUE)

```

## Versioning

#### Semantic versioning

`cytobox` uses a semantic versioning scheme, where package versions are of the form x.y.z:

- The patch version z is updated for very minor changes like bug fixes
- The minor version y is udpated for minor but important changes like addition of new functions
- The major version x is updated for major changes, which may introduce backwards incompatibiility

[![Build Status](https://travis-ci.com/fungenomics/cytobox.svg?branch=master)](https://travis-ci.com/fungenomics/cytobox)

# cytobox
A toolkit for analyzing single cell RNA-seq data created by the Kleinman Lab. The legacy package, `cytokit`,
which is no longer being developed, is available at https://github.com/sjessa/cytokit/. `cytobox` is
designed to work with [Seurat](https://satijalab.org/seurat/) objects and many functions are modeled on commonly used Seurat functions, to allow
for additional customization.

**Note**: `cytobox` currently works with Seurat V2.

## Installation

You will need to install [`devtools`](https://cran.r-project.org/web/packages/devtools/), and then run:

```r
devtools::install_github("fungenomics/cytobox", build_vignettes = TRUE)

```

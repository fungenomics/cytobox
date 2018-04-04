# Package-level documentation and package-level imports

#' cytokit: A toolkit for analyzing single cell RNA-seq data
#'
#' Shared Kleinman lab functions for processing, interrogating, and
#' visualizing single cell RNA-seq data.
#'
#' @docType package
#' @name cytokit
"_PACKAGE"

#' @importClassesFrom Seurat seurat
NULL

#' @importFrom magrittr %>%
NULL

#' pbmc
#'
#' @inherit Seurat::pbmc_small
"pbmc"

#' @import ggplot2
NULL

#' @importFrom ggplot2 theme
NULL

#' @importFrom glue glue
NULL

#' @importFrom viridis viridis
NULL

#' @import dplyr
NULL

#' @importFrom rlang :=
NULL

#' @importFrom stats median
NULL

#' @importFrom cowplot plot_grid
NULL

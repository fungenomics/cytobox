# Helper and utility functions

#' ggColours
#'
#' Get evenly spaced colours from around the colour wheel, which are the default
#' colours assigned to clusters by Seurat. The output of this function can be
#' passed to the \code{scale_colour_manual()} and \code{scale_fill_manual()} functions
#' from ggplot2, as the \code{values} argument. (\code{\link{ggColors}} points
#' to this function.)
#'
#' @param n Number of colours to return
#'
#' @return Named character vector, where names are the names of clusters, from
#' 0 to n-1, and values are the hex codes for the colours.
#' @export
#'
#' @examples
#'
#' n_clust <- 5
#' ggColours(n_clust)
#'
#' @references https://stackoverflow.com/a/8197703
#' @aliases ggColors
#' @importFrom grDevices hcl
ggColours <- function(n) {

    hues <- seq(15, 375, length = n + 1)
    colours <- hcl(h = hues, l = 65, c = 100)[1:n]
    names(colours) <- seq(0, n - 1) # Since the first cluster in Seurat is 0

    return(colours)

}

#' @export
ggColors <- ggColours


#' addEmbedding
#'
#' Given a Seurat object and a data frame where the rows correspond to cells,
#' in the same order as in the Seurat object, add two columns giving coordinates
#' in a dimensionality reduced space.
#'
#' @param object Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA() or Seurat::RunTSNE() to the object)
#' @param df Data frame where rows correspond to cells, in the same order
#' as in \code{rownames(object@@meta.data)}.
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves t-SNE by default. This should match the names of the elements of
#' the list object@@dr, so it will typically be one of "pca" or "tsne".
#' Default: "tsne"
#'
#' @export
#' @author Selin Jessa
#' @examples
#' df <- data.frame(Cell = rownames(pbmc@meta.data),
#'                  Cluster = pbmc@meta.data$res.0.8)
#'
#' addEmbedding(pbmc, df, reduction = "tsne")
addEmbedding <- function(object, df, reduction = "tsne") {

    # Get the axes for the reduced space
    # See here: http://dplyr.tidyverse.org/articles/programming.html#setting-variable-names
    vars <- colnames(object@dr[[reduction]]@cell.embeddings)[c(1, 2)]
    df %>%
        dplyr::mutate(!!vars[1] := object@dr[[reduction]]@cell.embeddings[, 1],
                      !!vars[2] := object@dr[[reduction]]@cell.embeddings[, 2])


}

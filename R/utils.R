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
#' @importFrom grDevices hcl hues
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
#' Given a Seurat seurat and a data frame where the rows correspond to cells,
#' in the same order as in the Seurat seurat, add two columns giving coordinates
#' in a dimensionality reduced space.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA() or Seurat::RunTSNE() to the object).
#' @param df Data frame with at least one column, "Cell", giving the cell ID
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves t-SNE by default. This should match the names of the elements of
#' the list seurat@@dr, so it will typically be one of "pca" or "tsne".
#' Default: "tsne"
#'
#' @export
#' @author Selin Jessa
#' @examples
#' df <- data.frame(Cell = rownames(pbmc@meta.data),
#'                  Cluster = pbmc@meta.data$res.1)
#'
#' addEmbedding(pbmc, df, reduction = "tsne")
addEmbedding <- function(seurat, df, reduction = "tsne") {

    # Get the axes for the reduced space
    # See here: http://dplyr.tidyverse.org/articles/programming.html#setting-variable-names
    vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(1, 2)]

    embedding <- data.frame(Cell = rownames(seurat@dr[[reduction]]@cell.embeddings)) %>%
        mutate(!!vars[1] := seurat@dr[[reduction]]@cell.embeddings[, 1],
               !!vars[2] := seurat@dr[[reduction]]@cell.embeddings[, 2])

    df <- dplyr::inner_join(df, embedding, by = "Cell")
    return(df)

}






#' fetchData
#'
#' Subset the seurat@@data matrix by gene and cluster.
#' Similar to Seurat::FetchData except it doesn't thrown an error if a gene
#' is not found in the data, and is more limited.
#'
#' @param seurat Seurat object
#' @param genes Genes to filter
#' @param clusters (Optional) Vector, include only cells with these identities
#' (e.g. cluster assignments). Searches in seurat@@ident.
#' @param return_cell Logical, whether or not to include a column with the cell ID.
#' Default: FALSE
#' @param return_cluster Logical, whether or not to include a column with the cluster.
#' Default: FALSE
#'
#' @return Expression matrix for genes specified
#' @export
#'
#' @author Selin Jessa
#' @examples
#' fetchData(pbmc, c("IL32", "MS4A1"))
#' fetchData(pbmc, c("IL32"), c(1, 2))
#' fetchData(pbmc, c("IL32", "MS4A1"), c(1, 2), return_cluster = TRUE, return_cell = TRUE)
fetchData <- function(seurat, genes, clusters = NULL, return_cell = FALSE, return_cluster = FALSE) {

    exp <- as.matrix(seurat@data)

    if(!any(genes %in% rownames(exp))) stop("No genes specified were ",
                                            "found in the data.")

    exp_filt <- as.data.frame(t(exp[which(rownames(exp) %in% genes),]))

    # Keep all
    if(is.null(clusters)) clusters <- unique(seurat@ident)
    ident_idx <- which(seurat@ident %in% clusters)

    # Handle only one gene case, and properly return a data frame
    if(nrow(exp_filt) == 1) {

        if(!is.null(ident)) exp_filt <- exp_filt[ident_idx]

        exp_filt <- as.data.frame(t(exp_filt))
        names(exp_filt) <- rownames(exp)[rownames(exp) %in% genes]

    } else {
        if(!is.null(ident)) exp_filt <- exp_filt[ident_idx,]
    }

    rownames(exp_filt) <- c() # Get rid of rownames

    if(return_cluster) exp_filt <- tibble::add_column(exp_filt, Cluster = seurat@ident[ident_idx], .before = 1)
    if(return_cell) exp_filt <- tibble::add_column(exp_filt, Cell = names(seurat@ident)[ident_idx], .before = 1)

    return(exp_filt)

}




#' theme_min
#'
#' A clean theme for ggplot2
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
#' @author Selin Jessa
#' @export
theme_min <- function(base_size = 11, base_family = "") {

    theme_light(base_size = 11, base_family = "") +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = "grey90", size = 1),
            strip.background = element_rect(fill = NA, colour = NA),
            strip.text.x = element_text(colour = "black", size = rel(1.2)),
            strip.text.y = element_text(colour = "black", size = rel(1.2)),
            title = element_text(size = rel(0.9)),
            axis.text = element_text(colour = "black", size = rel(0.8)),
            axis.title = element_text(colour = "black", size = rel(0.9)),
            legend.title = element_text(colour = "black", size = rel(0.9)),
            legend.key.size = unit(0.9, "lines"),
            legend.text = element_text(size = rel(0.7), colour = "black"),
            legend.key = element_rect(colour = NA, fill = NA),
            legend.background = element_rect(colour = NA, fill = NA)
        )
}


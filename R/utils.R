# Helper and utility functions





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
addEmbedding <- function(seurat, df, reduction = "tsne", dim1 = 1, dim2 = 2) {

    # Get the axes for the reduced space
    # See here: http://dplyr.tidyverse.org/articles/programming.html#setting-variable-names
    vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(dim1, dim2)]

    df$Cell <- as.character(df$Cell)

    embedding <- data.frame(Cell = seurat@cell.names, stringsAsFactors = FALSE) %>%
        mutate(!!vars[1] := seurat@dr[[reduction]]@cell.embeddings[, 1],
               !!vars[2] := seurat@dr[[reduction]]@cell.embeddings[, 2])

    df <- dplyr::inner_join(df, embedding, by = "Cell")
    return(df)

}




#' findGenes
#'
#' Given a set of genes, find the ones which are detected in the sample,
#' and which are not.
#'
#' @param seurat Seurat object
#' @param genes Character vector of genes
#'
#' @return A named list with two elements: "detected" and "undetected"
#' each storing character vectors with the genes in each category
#' @export
#'
#' @author Selin Jessa
#' @examples
#' find_out <- findGenes(pbmc, c("IL32", "CD79B", "foo"))
#' find_out$detected
#' find_out$undetected
findGenes <- function(seurat, genes) {

    genes_detected <- genes[genes %in% rownames(seurat@data)]
    genes_undetected <- setdiff(genes, genes_detected)

    list(detected = genes_detected, undetected = genes_undetected)

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
#' @param scaled Logical, whether to fetch scaled data from seurat@@scale.data.
#' Default: FALSE.
#'
#' @return Expression matrix for genes specified
#' @export
#'
#' @author Selin Jessa
#' @examples
#' fetchData(pbmc, c("IL32", "MS4A1"))
#' fetchData(pbmc, c("IL32"), c(1, 2))
#' fetchData(pbmc, c("IL32", "MS4A1"), c(1, 2), return_cluster = TRUE, return_cell = TRUE)
fetchData <- function(seurat, genes, clusters = NULL,
                      return_cell = FALSE, return_cluster = FALSE, scaled = FALSE) {

    genes_out <- findGenes(seurat, genes)

    n_undetected <- length(genes_out$undetected)

    if (n_undetected > 0) {

        if (n_undetected > 10) {

            message(paste0("NOTE: [",
                         paste0(head(genes_out$undetected, 10), collapse = ", "),
                         "] and ", n_undetected - 10, " other genes are undetected in ", seurat@project.name))

        } else {

            message(paste0("NOTE: [",
                         paste0(genes_out$undetected, collapse = ", "),
                         "] undetected in ", seurat@project.name))

        }
    }

    if (length(genes_out$detected) == 0) stop("No genes specified were ",
                                             "found in the data.")

    if (scaled) exp <- as.matrix(seurat@scale.data)
    else exp <- as.matrix(seurat@data)

    exp_filt <- as.data.frame(t(exp[which(rownames(exp) %in% genes_out$detected),]))

    # Keep all
    if(is.null(clusters)) clusters <- levels(seurat@ident)
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




#' Retrieve cells in specified clusters
#'
#' Given a Seurat object and the names of specific clusters (corresponding
#' to identities in \code{seurat@@ident}), return the names of cells
#' in those clusters.
#'
#' @param seurat Seurat object
#' @param clusters Vector specifying clusters for which cells should be retrieved.
#' Elements should correspond to levels in \code{seurat@@ident}.
#'
#' @return Character vector of cell names
#' @export
#' @author Selin Jessa
#'
#' @examples
#' whichCells(pbmc, clusters = c(0, 1))
whichCells <- function(seurat, clusters) {

    names(seurat@ident)[seurat@ident %in% clusters]

}




#' getVarianceExplained
#'
#' Compute variance explained by PCA, given a Seurat object for which the PCA
#' dim. reduction has been calculated
#'
#' @param seurat Seurat object
#' @param n Numeric, number of PCs for which variance should be reported.
#' Default: 10
#'
#' @return List with two vectors, "percent.var.explained" and "cumulative.var.explained",
#' reported for the first \code{n} PCs
#'
#' @examples
#' getVarianceExplained(pbmc, n = 5)
#'
#' @author Adapted from Yang Yang
#' @export
getVarianceExplained <- function(seurat, n = 10) {

    sdev <- seurat@dr$pca@sdev
    variance <- sdev^2
    sum.variance <- sum(variance)
    proportion.variance <- variance/sum.variance * 100
    acc_prop_var <- cumsum(proportion.variance)

    return(list(percent.var.explained = head(proportion.variance, n),
                cum.var.explained = head(acc_prop_var, n)))

}

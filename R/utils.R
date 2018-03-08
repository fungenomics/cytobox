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
addEmbedding <- function(seurat, df, reduction = "tsne") {

    # Get the axes for the reduced space
    # See here: http://dplyr.tidyverse.org/articles/programming.html#setting-variable-names
    vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(1, 2)]

    df$Cell <- as.character(df$Cell)

    embedding <- data.frame(Cell = seurat@cell.names, stringsAsFactors = FALSE) %>%
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
    if (length(genes_out$undetected > 0)) message(paste0("NOTE: [",
                                                         paste0(genes_out$undetected, collapse = ", "),
                                                         "] undetected in the data"))

    if(length(genes_out$detected) == 0) stop("No genes specified were ",
                                             "found in the data.")

    if (scaled) exp <- as.matrix(seurat@scale.data)
    else exp <- as.matrix(seurat@data)

    exp_filt <- as.data.frame(t(exp[which(rownames(exp) %in% genes_out$detected),]))

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






# Functions for studying correlation of gene expression



#' meanClusterExpression
#'
#' This function computes the mean expression of each gene in each cluster
#'
#' @param seurat Seurat object
#' @param genes (Optional) Character vector specifying genes for which means
#' should be calculated. If any provided genes are not detected, they are
#' silently dropped. Default: all detected genes, i.e. \code{rownames(seurat@@data)}.
#'
#' @return A gene x cluster matrix, where each column contains the mean
#' expression of each gene over all the cells in that cluster
#' @export
#' @author Selin Jessa
#'
#' @examples
#' meanClusterExpression(pbmc, c("IL32", "CD2"))
meanClusterExpression <- function(seurat, genes = NULL) {

    meanClusterExpression.percluster <- function(cluster, seurat, genes) {

        cells <- cytokit::whichCells(seurat, cluster)
        expr <- as.matrix(seurat@data[genes, cells])
        rowMeans(expr)

    }

    message(glue("Computing cluster means for {seurat@project.name}..."))

    if (is.null(genes)) genes <- rownames(seurat@data)
    genes <- findGenes(seurat, genes) %>% .$detected

    expr_means <- do.call(rbind, lapply(levels(seurat@ident),
                                        meanClusterExpression.percluster,
                                        seurat,
                                        genes)) %>%
        magrittr::set_rownames(levels(seurat@ident)) %>%
        t()

    return(expr_means)

}



#' Correlate expression
#'
#' Compute a correlation matrix of gene expression. For matrices, this function will
#' correlate the columns; for Seurat objects, this function automatically retrieves
#' the cluster means of the provided genes using \code{\link{meanClusterExpression}},
#' and then correlates the cluster means.
#'
#' @param s1 A Seurat object, or a gene x cluster/sample matrix giving mean
#' cluster/sample expression for dataset 1. If a matrix, rownames should be gene symbols. This dataset will form the x-axis of the plot;
#' column names will be pulled from the column names of this object.
#' @param s2 A Seurat object, or a gene x cluster/sample matrix giving mean cluster/sample
#' expression for dataset 2.
#' This dataset will form the y-axis of the plot; row names will be pulled from the
#' column names of this object.
#' @param genes Character vector of genes, in the same species as \code{s1}.
#' @param from_sp Species of \code{s1}
#' @param to_sp Species of \code{s2}. If different from \code{from_sp}, conversion
#' of gene names will be handled automatically.
#' @param method Character, one of "pearson", "kendall", "spearman", specifying
#' the method to use for the correlations. Default: "spearman"
#'
#' @export
#' @author Selin Jessa
#'
#'
#' @examples
#' # Compute pairwise correlations between clusters in pbmc
#' correlateExpression(s1 = pbmc,
#'                     s2 = pbmc,
#'                     genes = head(rownames(pbmc@@raw.data), 100),
#'                     from_sp = "hg",
#'                     to_sp = "hg")
correlateExpression <- function(s1, s2, genes, from_sp, to_sp,
                                return_input = FALSE) {

    if (class(s1) == "seurat") mat1 <- meanClusterExpression(s1, genes)
    else mat1 <- s1

    if (class(s2) == "seurat") {

        # If species are not the same, don't use the provided genes for s2
        if (from_sp != to_sp) mat2 <- meanClusterExpression(s2)
        else mat2 <- meanClusterExpression(s2, genes)

    } else mat2 <- s2

    # 1. If needed, convert matrices to have gene names in the same species (as s1)
    if (from_sp != to_sp) {

        # mat2 is in hg, convert to mm
        if (to_sp == "hg") mat2 <- hg2mm(mat2)
        # mat2 is in mm, convert to hg
        else if (to_sp == "mm") mat2 <- mm2hg(mat2)

    }
    # Now we should have two matrices with gene names in the same species

    # 2. Subset countmats to genes present in both datasets
    # (i.e. accounting for ones that are not detected in one sample or the other,
    # and genes that may have gotten lost in translation between species)
    genes_keep <- genes[(genes %in% rownames(mat1)) & (genes %in% rownames(mat2))]

    mat1_sub <- mat1[genes_keep, ]
    mat2_sub <- mat2[genes_keep, ]

    if (return_input) return(list(mat1_sub, mat2_sub))

    # 3. Check
    if (all(rownames(mat1_sub) == rownames(mat2_sub))) {

        message("Computing correlations...")

        # 4. Take the correlation
        cormat <- cor(mat1_sub, mat2_sub, method = "spearman", use = "complete.obs")
        return(cormat)

    } else {

        stop("Error in subsetting count matrices, rownames are not equal")

    }

}

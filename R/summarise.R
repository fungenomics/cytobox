# Functions for computing summary stats


#' meanMarkerExpression
#'
#' @param seurat Seurat object
#' @param markers String or character vector specifying gene(s) to use
#'
#' @return Data frame with two columns: "Cell" which specifies the cell ID
#' obtained from \code{colnames(seurat@@data)}, and "Mean_marker_expression"
#' which is the expression value for the marker gene, or mean if multiple
#' genes were provided.
#'
#' @export
#' @author Selin Jessa
meanMarkerExpression <- function(seurat, markers) {

    # Get expression data from the Seurat seurat
    exp <- as.data.frame(as.matrix(seurat@data))

    # Filter to the markers
    filt_exp <- exp[rownames(exp) %in% markers, ]

    # Calculate mean expression
    mean_exp <- colSums(filt_exp)/nrow(filt_exp)
    exp_df <- data.frame(Cell = colnames(exp),
                         Mean_marker_expression = as.numeric(mean_exp))

    return(exp_df)
}


#' perecentilesMarkerExpression
#'
#' @param seurat Seurat object
#' @param markers String or character vector specifying gene(s) to use
#'
#' @return Data frame with three columns: "Cell" which specifies the cell ID
#' obtained from \code{colnames(seurat@@data)}, "Percentile" giving the
#' percentile expression for each cell, and
#' "Gradient_group" which gives the ordering that should be used to assign
#' colours by gradient.
#'
#' @export
#' @author adapted from code from Alexis Blanchet-Cohen
#' @importFrom stats ecdf
percentilesMarkerExpression <- function(seurat, markers) {

    exp.matrix <- as.data.frame(as.matrix(seurat@data)) %>%
        tibble::rownames_to_column(var = "gene_name")

    exp.matrix.filt <- filter(exp.matrix, gene_name %in% markers)
    exp.matrix.filt.summed <- colSums(select(exp.matrix.filt, -gene_name))
    # Notice that cell.type is really the sum of expression values for all the genes
    # in the markers argument
    exp.matrix.filt.summed <- data.frame(cell.type = unname(exp.matrix.filt.summed),
                                          Cell = colnames(seurat@data))

    exp <- exp.matrix.filt.summed

    # Add percentiles
    exp.zero     <- filter(exp, cell.type == 0) %>% mutate(percentile=NA)
    exp.not.zero <- filter(exp, cell.type != 0)
    exp.not.zero$percentile <- ecdf(exp.not.zero$cell.type)(exp.not.zero$cell.type) * 100
    exp <- rbind(exp.zero, exp.not.zero)

    # Add color columns
    exp <- exp %>%
        dplyr::mutate(Gradient_group = case_when(
            percentile > 0  & percentile <= 50  ~ 2,
            percentile > 50 & percentile <= 70  ~ 3,
            percentile > 70 & percentile <= 90  ~ 4,
            percentile > 90 & percentile <= 92  ~ 5,
            percentile > 92 & percentile <= 94  ~ 6,
            percentile > 94 & percentile <= 96  ~ 7,
            percentile > 96 & percentile <= 98  ~ 8,
            percentile > 98 & percentile <= 100 ~ 9,
            TRUE ~ 1)) %>%
        dplyr::select(Cell, Cell.type = cell.type, Percentile = percentile,
                      Gradient_group)

    exp <- dplyr::arrange(exp, Cell.type)

    return(exp)

}

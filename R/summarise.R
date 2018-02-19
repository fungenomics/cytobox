# Functions for computing summary stats


#' meanMarkerExpression
#'
#' @param object Seurat object
#' @param markers String or character vector specifying gene(s) to use
#'
#' @return Data frame with two columns: "Cell" which specifies the cell ID
#' obtained from \code{colnames(object@@data)}, and "Mean_marker_expression"
#' which is the expression value for the marker gene, or mean if multiple
#' genes were provided.
#'
#' @export
#' @author Selin Jessa
meanMarkerExpression <- function(object, markers) {

    # Get expression data from the Seurat object
    exp <- object@data %>% as.matrix %>% as.data.frame

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
#' @param object Seurat object
#' @param markers String or character vector specifying gene(s) to use
#'
#' @return Data frame with three columns: "Cell" which specifies the cell ID
#' obtained from \code{colnames(object@@data)}, "Percentile" giving the
#' percentile expression for each cell, and
#' "Gradient_group" which gives the ordering that should be used to assign
#' colours by gradient.
#'
#' @export
#' @author adapted from code from Alexis Blanchet-Cohen
#' @importFrom stats ecdf
percentilesMarkerExpression <- function(object, markers) {

    # Get expression data
    exp <- object@data %>%
        as.matrix() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "gene_name")

    exp_filt <- exp %>% filter(gene_name %in% markers)
    exp_filt_summed <- exp_filt %>% dplyr::select(-gene_name) %>% colSums()
    exp_filt_summed <- data.frame(Cell = colnames(object@data),
                                  Exp_sum = unname(exp_filt_summed))

    exp_zero <- exp_filt_summed %>% filter(Exp_sum == 0) %>%
        dplyr::mutate(percentile = NA)
    exp_not_zero <- exp_filt_summed %>% filter(Exp_sum != 0)

    # Add percentiles
    exp_not_zero$percentile <- ecdf(exp_not_zero$Exp_sum)(exp_not_zero$Exp_sum) * 100
    exp_data <- rbind(exp_zero, exp_not_zero)

    # Add colour columns based on percentiles
    exp_data <- exp_data %>%
        dplyr::mutate(gradient_group = case_when(
            percentile > 0  & percentile <= 50  ~ 2,
            percentile > 50 & percentile <= 70  ~ 3,
            percentile > 70 & percentile <= 90  ~ 4,
            percentile > 90 & percentile <= 92  ~ 5,
            percentile > 92 & percentile <= 94  ~ 6,
            percentile > 94 & percentile <= 96  ~ 7,
            percentile > 96 & percentile <= 98  ~ 8,
            percentile > 98 & percentile <= 100 ~ 9,
            TRUE ~ 1
        )) %>%
        dplyr::arrange(Exp_sum) %>%
        dplyr::select(Cell, Percentile = percentile, Gradient_group = gradient_group)

    return(exp_data)

}

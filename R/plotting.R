# Functions for plotting


#' tsneByMeanMarkerExpression
#'
#' Plot a low-dimensional embedding of the cells,
#' coloured by expression of a gene, or mean expression of a group of marker
#' genes. Defaults to t-SNE space, but see the \code{reduction} argument for
#' how to plot in PCA space instead.
#'
#' @param object Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA() or Seurat::RunTSNE() to the object)
#' @param markers String or character vector specifying gene(s) to use
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves t-SNE by default. This should match the names of the elements of
#' the list object@@dr, so it will typically be one of "pca" or "tsne".
#' Default: "tsne"
#'
#' @export
#'
#' @return A ggplot object
#'
#' @author Selin Jessa
#' @examples
#' tsneByMeanMarkerExpression(pbmc, "IL32")
#' tsneByMeanMarkerExpression(pbmc, c("IL32", "CD2"), reduction = "pca")
tsneByMeanMarkerExpression <- function(object, markers,
                                       reduction = "tsne") {

    # Get mean expression for markers
    exp_df <- meanMarkerExpression(object, markers)

    # Get dimensionality reduction coordinates
    exp_df <- object %>% addEmbedding(exp_df, reduction)

    # Get the variable names
    vars <- colnames(object@dr[[reduction]]@cell.embeddings)[c(1, 2)]

    # Plot
    gg <- exp_df %>%
        dplyr::arrange(Mean_marker_expression) %>% # Order in which points will be plot
        ggplot(aes(x = exp_df[[vars[1]]], y = exp_df[[vars[2]]])) +
        geom_point(aes(colour = Mean_marker_expression), size = rel(0.8)) +
        viridis::scale_color_viridis() +
        xlab(vars[1]) + ylab(vars[2])

    gg <- gg +
        theme_min()

    return(gg)

}


#' tsneByPercentileMarkerExpression
#'
#' Plot a low-dimensional embedding of the cells,
#' coloured by expression percentile of a gene, of the total expression of a
#' group of marker genes.
#'
#' @param object Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA() or Seurat::RunTSNE() to the object)
#' @param markers String or character vector specifying gene(s) to use
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves t-SNE by default. This should match the names of the elements of
#' the list object@dr, so it will typically be one of "pca" or "tsne".
#' Default: "tsne"
#'
#' @export
#'
#' @return A ggplot object
#'
#' @author Selin Jessa
#' @examples
#' tsneByPercentileMarkerExpression(pbmc, "IL32")
#' tsneByPercentileMarkerExpression(pbmc, c("IL32", "CD2"), reduction = "pca")
tsneByPercentileMarkerExpression <- function(object, markers, reduction = "tsne") {

    # Get expression percentiles
    percentiles <- percentilesMarkerExpression(object, markers)

    color_grad_labels <- c("Undetected",
                           "> 0 & \u2265 50",
                           "> 50 & \u2265 70",
                          "> 70 & \u2265 90",
                          "> 90 & \u2265 92",
                          "> 92 & \u2265 94",
                          "> 94 & \u2265 96",
                          "> 96 & \u2265 98",
                          "> 98 & \u2265 100")

    # Get the variable names, and embedding
    vars <- colnames(object@dr[[reduction]]@cell.embeddings)[c(1, 2)]
    percentiles <- object %>% addEmbedding(percentiles, reduction)

    # Plot
    gg <- percentiles %>%
        dplyr::arrange(Percentile) %>%
        ggplot(aes(x = percentiles[[vars[1]]], y = percentiles[[vars[2]]])) +
        geom_point(aes(colour = factor(Gradient_group, levels = seq(1, 9))),
                   alpha = 0.6) +
        scale_color_manual(values = viridis(9),
                           name = "Expression level percentile",
                           labels = color_grad_labels,
                           # Ensured all levels are displayed in the legend
                           drop = FALSE) +
        guides(colour = guide_legend(order = 1)) +
        xlab(vars[1]) + ylab(vars[2])

    gg <- gg + theme_min()

    return(gg)

}



#' @importFrom ggplot2 theme_light
NULL


#' theme_min
#'
#' A clean theme for ggplot2
#'
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

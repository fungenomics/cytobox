# Functions for general plotting


#' tsneByMeanMarkerExpression
#'
#' Plot a low-dimensional embedding of the cells,
#' coloured by expression of a gene, or mean expression of a group of marker
#' genes. Defaults to t-SNE space, but see the \code{reduction} argument for
#' how to plot in PCA space instead.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA() or Seurat::RunTSNE() to the object)
#' @param markers String or character vector specifying gene(s) to use
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves t-SNE by default. This should match the names of the elements of
#' the list seurat@@dr, so it will typically be one of "pca" or "tsne".
#' Default: "tsne"
#'
# @export
#'
#' @return A ggplot object
#'
#' @author Selin Jessa
# @examples
# tsneByMeanMarkerExpression(pbmc, "IL32")
# tsneByMeanMarkerExpression(pbmc, c("IL32", "CD2"), reduction = "pca")
tsneByMeanMarkerExpression <- function(seurat, markers,
                                       reduction = "tsne") {

    # Get mean expression for markers
    exp_df <- meanMarkerExpression(seurat, markers)

    # Get dimensionality reduction coordinates
    exp_df <- seurat %>% addEmbedding(exp_df, reduction)

    # Get the variable names
    vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(1, 2)]

    # Plot
    gg <- exp_df %>%
        dplyr::arrange(Mean_marker_expression) %>% # Order in which points will be plot
        ggplot(aes(x = exp_df[[vars[1]]], y = exp_df[[vars[2]]])) +
        geom_point(aes(colour = Mean_marker_expression), size = rel(0.8), alpha = 0.6) +
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
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA() or Seurat::RunTSNE() to the object)
#' @param markers String or character vector specifying gene(s) to use
#' @param label Logical, whether to label clusters on the plot. Default: TRUE.
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves t-SNE by default. This should match the names of the elements of
#' the list seurat@@dr, so it will typically be one of "pca" or "tsne".
#' Default: "tsne"
#' @param palette String, one of "viridis" or "blues", specifying which gradient
#' palette to use. Default: viridis.
#'
#' @export
#'
#' @return A ggplot object
#'
#' @author Selin Jessa
#' @examples
#' tsneByPercentileMarkerExpression(pbmc, "IL32")
#' tsneByPercentileMarkerExpression(pbmc, c("IL32", "CD2"), reduction = "pca")
tsneByPercentileMarkerExpression <- function(seurat, markers,
                                             label = TRUE,
                                             reduction = "tsne",
                                             palette = "viridis",
                                             extra = FALSE) {

    # Get expression percentiles
    percentiles <- percentilesMarkerExpression(seurat, markers)

    color_grad_labels <- c("Undetected",
                           "> 0 & \u2264 50",
                           "> 50 & \u2264 70",
                          "> 70 & \u2264 90",
                          "> 90 & \u2264 92",
                          "> 92 & \u2264 94",
                          "> 94 & \u2264 96",
                          "> 96 & \u2264 98",
                          "> 98 & \u2264 100")
    alpha_grad <- c(0.025, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
    alpha_grad_labels <- color_grad_labels

    # TODO you can just set the alpha group to the gradient group,
    # and as its values, pass the discrete gradient by hand... saves some code

    # Get the variable names, and embedding
    vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(1, 2)]
    percentiles <- addEmbedding(seurat, percentiles, reduction)

    # Plot
    gg <- ggplot(percentiles,
                 aes(x = percentiles[[vars[1]]], y = percentiles[[vars[2]]])) +
        geom_point(aes(colour = factor(Gradient_group, levels = seq(1, 9)),
                       alpha = factor(Gradient_group, levels = seq(1, 9))))

    if (palette == "viridis") {

        gg <- gg +
            scale_color_manual(
                values = viridis(9),
                name = "Expression level percentile",
                labels = color_grad_labels,
                # Ensure all levels are displayed in the legend
                drop = FALSE)

    } else if (palette == "blues") {

        gg <- gg +
            scale_color_manual(
                values = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
                name = "Expression level percentile",
                labels = color_grad_labels,
                drop = FALSE)
    }

    gg <- gg +
        guides(colour = guide_legend(order = 1)) +
        scale_alpha_manual(values = alpha_grad, labels = alpha_grad_labels,
                           name = "Expression level percentile", drop = FALSE) +
        xlab(vars[1]) + ylab(vars[2]) +
        theme_min()

    if (label) {

        if (reduction == "pca") {

            message("Plotting labels is currently only available for reduction = 'tsne';",
                    " returning plot without labels.")
            return(gg)
        }

        centers <- clusterCenters(seurat)

        gg <- gg +
            ggrepel::geom_label_repel(data = centers,
                             aes(x = mean_tSNE_1, y = mean_tSNE_2),
                             label = centers$Cluster,
                             segment.color = 'grey50',
                             fontface = 'bold',
                             alpha = 0.8,
                             segment.alpha = 0.8,
                             label.size = NA,
                             force = 2,
                             nudge_x = 5, nudge_y = 5,
                             segment.size = 0.5,
                             arrow = arrow(length = unit(0.01, 'npc')))

    }

    return(gg)

}


#' theme_min
#'
#' A clean theme for ggplot2
#'
#' @importFrom ggplot2 theme_light
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

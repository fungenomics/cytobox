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
#' @param title (Optional) string used as title for the plot.
#' @param palette String, one of "viridis" or "blues", specifying which gradient
#' palette to use. Default: viridis.
#' @param extra Logical, plot a detailed three-panel plot, where the first
#' is a proportional bar plot of cells in each cluster in each percentile
#' group, the second is a ridge plot showing density in each cluster of the mean
#' expression of the markers (coloured by median percentile group within the
#' cluster), and the third is the labelled tSNE plot coloured
#' by percentiles. Requires \code{label = TRUE}. Default: FALSE.
#' Default: FALSE.
#' @param verbose Logical, whether to print status updates. Default: FALSE.
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
                                             title = NULL,
                                             palette = "viridis",
                                             extra = FALSE,
                                             verbose = FALSE) {

    if (verbose) message("Computing percentiles...")

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

    if (verbose) message("Constructing percentiles plot...")

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

        if (extra) {

            # Frequency plot
            clusters <- data.frame(Cell = names(seurat@ident), Cluster = seurat@ident)
            df <- full_join(percentiles, clusters, by = "Cell")

            if (verbose) message("Computing frequencies...")

            freq <- df %>%
                mutate(Gradient_group = factor(Gradient_group, levels = seq(9, 1)),
                       Cluster = factor(Cluster)) %>%
                group_by(Cluster, Gradient_group) %>%
                summarise(n = n()) %>%
                group_by(Cluster) %>%
                mutate(Proportion = n/sum(n)) %>%
                ungroup() %>%
                group_by(Cluster) %>%
                mutate(n_per_clust = sum(n)) %>%
                ungroup()

            # Get the order of clusters, ranked by proportion of cells
            # in the top percentile group
            rank <- freq %>%
                arrange(Gradient_group, desc(Proportion)) %>%
                group_by(Cluster) %>%
                slice(1) %>%
                arrange(desc(Proportion)) %>%
                ungroup() %>%
                mutate(rank = seq(1:length(unique(freq$Cluster)))) %>%
                select(Cluster, rank)

            # Number of cells in each cluster
            n_cells <- freq %>% group_by(Cluster) %>% summarise(n_cells_in_clust = sum(n))

            freq2 <- freq %>%
                left_join(n_cells, by = "Cluster") %>%
                left_join(rank, by = "Cluster") %>%
                arrange(rank) %>%
                mutate(Cluster = glue("{Cluster} (n={n_cells_in_clust})"))

            freq2$Cluster <- factor(freq2$Cluster,
                                    levels = rev(unique(freq2$Cluster)))

            if (verbose) message("Constructing bar plot...")

            p1 <- freq2 %>%
                ggplot(aes(x = Cluster, y = Proportion)) +
                geom_bar(aes(fill = Gradient_group), stat = "identity", position = "stack",
                         width = 0.4) +
                scale_fill_manual(values = rev(viridis(9)),
                                  name = "Expression level percentile",
                                  labels = rev(color_grad_labels),
                                  # Ensure all levels are displayed in the legend
                                  drop = FALSE) +
                theme_min() +
                coord_flip() +
                theme(legend.position = "none")

            if (verbose) message("Computing means...")

            # Ridge plot
            df_means <- meanMarkerExpression(seurat, markers) %>%
                full_join(df, by = "Cell") %>%
                group_by(Cluster) %>%
                mutate(Median = median(Gradient_group)) %>%
                mutate(Median = factor(Median, levels = seq(9, 1)))

            if (verbose) message("Constructing ridge plot...")

            p2 <- df_means %>% ggplot(aes(x = Mean_marker_expression,
                                          y = factor(Cluster, levels = rev(rank$Cluster)))) +
                ggridges::geom_density_ridges(scale = 0.7, rel_min_height = 0.0005, size = 0.3,
                                              aes(fill = Median)) +
                scale_fill_manual(values = rev(viridis(9)),
                                  name = "Expression level percentile",
                                  labels = rev(color_grad_labels),
                                  # Ensure all levels are displayed in the legend
                                  drop = FALSE) +
                theme_min() +
                scale_x_continuous(expand = c(0.01, 0)) +
                scale_y_discrete(expand = c(0.01, 0)) +
                xlab("Mean expression") +
                theme(legend.position = "none",
                      axis.title.y = element_blank(),
                      axis.text.y  = element_blank(),
                      axis.ticks.y = element_blank())

            if (verbose) message("Combining plots...")

            if (!is.null(title)) {

                p1 <- p1 + ggtitle(title)
                p2 <- p2 + ggtitle("") # So that all plots are aligned
                gg <- gg + ggtitle("")

            }

            combined <- cowplot::plot_grid(p1, p2, gg, rel_widths = c(0.3, 0.2, 1), nrow = 1)
            return(combined)

        }

    }

    if (!is.null(title)) {

        gg <- gg + ggtitle(title)
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

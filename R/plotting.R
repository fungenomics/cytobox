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
#' @param genes String or character vector specifying gene(s) to use
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
tsneByMeanMarkerExpression <- function(seurat, genes,
                                       reduction = "tsne") {

    # Get mean expression for markers
    exp_df <- meanMarkerExpression(seurat, genes)

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
#' @param genes String or character vector specifying gene(s) to use
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
#' @aliases dashboard
#' @examples
#' tsneByPercentileMarkerExpression(pbmc, "IL32")
#' tsneByPercentileMarkerExpression(pbmc, c("IL32", "CD2"), reduction = "pca")
tsneByPercentileMarkerExpression <- function(seurat, genes,
                                             label = TRUE,
                                             reduction = "tsne",
                                             title = NULL,
                                             palette = "viridis",
                                             extra = FALSE,
                                             verbose = FALSE) {

    if (verbose) message("Computing percentiles...")

    # Get expression percentiles
    percentiles <- percentilesMarkerExpression(seurat, genes)

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
            df_means <- meanMarkerExpression(seurat, genes) %>%
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


#' @export
dashboard <- function(seurat, genes,
                      title = NULL,
                      verbose = FALSE) {


    tsneByPercentileMarkerExpression(seurat, genes, title = title,
                                     extra = TRUE, verbose = verbose)

}




#' vln
#'
#' Similar to Seurat::VlnPlot() except it silently skips genes that are not found in the
#' data, and has an option to group plots by cluster instead of gene.
#'
#' @param seurat Seurat object
#' @param genes Genes to plot violins for
#' @param facet_by String, one of "gene" or "cluster". Default: "cluster",
#' genes will be on the x axis, and there will be one plot per cluster.
#' If "gene", clusters will be on the x axis, and there will be one plot per
#' gene (akin to Seurat::VlnPlot)
#'
#' @return A ggplot2 object
#' @author Selin Jessa
#' @export
#'
#' @examples
#' vln(pbmc, c("IL32", "MS4A1"))
#' vln(pbmc, c("IL32", "MS4A1"), facet_by = "gene")
vln <- function(seurat, genes, facet_by = "cluster") {

    exp <- fetchData(seurat, genes, return_cluster = TRUE) %>%
        tidyr::gather(Gene, Expression, 2:length(.))

    if (facet_by == "gene") gg <- ggplot(exp, aes(x = Cluster, y = Expression, fill = Cluster))
    else if (facet_by == "cluster") gg <- ggplot(exp, aes(x = Gene, y = Expression, fill = Cluster))

    gg <- gg +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        geom_jitter(size = 0.4, alpha = 0.5) +
        theme_min() +
        theme(legend.position = "none")

    if (facet_by == "gene") {

        gg <- gg + facet_wrap(~ Gene, scales = "free_y", ncol = 2) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

    } else if (facet_by == "cluster") gg <- gg + facet_wrap(~ Cluster, scales = "free_y")

    return(gg)

}


#' vlnGrid
#'
#' A grid of small violin plots, one plot per gene per cluster, similar to
#' Figure 5D in the Drop-seq paper by Macosko et al.
#'
#' @param seurat Seurat object
#' @param genes Genes to plot violins for
#'
#' @return A ggplot2 object
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
# Use the first 15 genes in the data
#' vlnGrid(pbmc, head(rownames(pbmc@data), 15))
vlnGrid <- function(seurat, genes) {

    expr <- fetchData(seurat, genes, return_cell = TRUE, return_cluster = TRUE) %>%
        tidyr::gather(Marker, Expression, 3:length(.))

    expr$Cluster <- factor(expr$Cluster, levels = seq(length(unique(seurat@ident))-1, 0))

    gg <- expr %>%
        ggplot(aes(x = Cluster, y = Expression)) +
        geom_violin(aes(fill = Cluster), scale = "width", size = 0.5) +
        facet_wrap(~ Marker, ncol = length(unique(expr$Marker))) +
        theme_min() +
        coord_flip() +
        ggplot2::theme(panel.border = element_blank(),
              axis.ticks.x = element_blank(),
              axis.text.x = element_blank(),
              axis.line.x = element_blank(),
              legend.position = "none",
              strip.text.x = element_text(angle = 30, size = 8))

    return(gg)

}




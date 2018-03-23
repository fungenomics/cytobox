# Functions for general plotting




#' tsne
#'
#' Plot a coloured, labelled tSNE plot for a datasets, akin to Seurat::TSNEPlot()
#'
#' @param seurat Seurat object, where Seurat::RunTSNE() has been applied
#' @param colour_by (Optional) String, specifying the column in \code{seurat@@meta.data}
#' by which to colour cells. Default: NULL, colour cells by cluster (in \code{seurat@@ident}).
#' @param colours (Optional) Character vector of colours for points. If \code{colour_by}
#' is NULL, cells will be coloured by cluster; this should be a named character vector of colours for points. Names should
#' correspond to cluster names (e.g. \code{levels(seurat@@ident)}). If
#' specifying \code{colour_by}, and the variable is discrete, this should be a named vector corresponding to categorical
#' values in the column of \code{seurat@@meta.data} specified. Otherwise, if the variable
#' is continuous, pass the gradient to use, or a few colours (from low to high) from which a gradient
#' should be created, and specify \code{colour_by_type = "continuous"}. The default is to
#' use ggplot2 colours.
#' @param colour_by_type (Optional) String, one of "discrete" or "continuous".
#' If specifying \code{colour_by} and providing colours to the \code{colours}
#' argument, specify whether the \code{colour_by} variable is discrete or continuous.
#' Default: discrete.
#' @param label Logical, whether to plot cluster labels. Default: TRUE
#' @param point_size Numeric, size of points in scatter plot. Default: 0.6
#' @param alpha Numeric, fixed alpha value for points: Default: 0.8
#' @param legend Logical, whether to plot legend. Default: FALSE if \code{colour_by}
#' is NULL and \code{label} is TRUE, true otherwise.
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param title (Optional) String specifying title.
#' @param hide_ticks Logical, whether to hide axis ticks. Default: FALSE
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at seurat@@ident) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#'
#' @return A ggplot2 object
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' tsne(pbmc)
tsne <- function(seurat,
                 colour_by = NULL,
                 colours = NULL,
                 colour_by_type = "discrete",
                 label = TRUE, point_size = 0.6, alpha = 0.8,
                 legend = ifelse((is.null(colour_by)) && (label), FALSE, TRUE),
                 label_repel = TRUE,
                 label_size = 4,
                 hide_ticks = FALSE,
                 title = NULL,
                 label_short = FALSE) {

    embedding <- data.frame(Cell = seurat@cell.names,
                            tSNE_1 = seurat@dr$tsne@cell.embeddings[, 1],
                            tSNE_2 = seurat@dr$tsne@cell.embeddings[, 2],
                            stringsAsFactors = FALSE)

    if (label && all(is.na(seurat@ident))) {
        label <- FALSE
        message("NOTE: identity of all cells is NA, setting 'label' to FALSE.")
    }

    if (is.null(colour_by)) {

        embedding$Cluster <- seurat@ident
        gg <- ggplot(embedding, aes(x = tSNE_1, y = tSNE_2))

        if (is.null(colours)) {

            # Assuming that the order of the levels is correct in the seurat object,
            # this should find the colours of the original clusters, and whatever they've been renamed,
            # if and only if the number of new cluster IDs is equal to the number of old ones
            colours <- ggColors(length(levels(seurat@ident)))
            names(colours) <- levels(seurat@ident)

        }

        gg <- gg +
            geom_point(aes(colour = Cluster), size = point_size, alpha = alpha) +
            scale_color_manual(values = colours)

    } else {

        embedding[[colour_by]] <- seurat@meta.data[[colour_by]]
        gg <- ggplot(embedding, aes(x = tSNE_1, y = tSNE_2))

        gg <- gg +
            geom_point(aes_string(colour = colour_by), size = point_size, alpha = alpha)


        if (!is.null(colours)) { # Otherwise default ggplot2 colours are used

            if (colour_by_type == "discrete") gg <- gg + scale_color_manual(values = colours)
            else if (colour_by_type == "continuous") gg <- gg + scale_color_gradientn(colours = colours)
        }

    }

    if (label) {

        centers <- clusterCenters(seurat)
        gg <- gg + addLabels(centers, label_repel, label_size, label_short)

    }

    gg <- gg + theme_min() + xlab("tSNE 1") + ylab("tSNE 2")

    if (!legend) gg <- gg + noLegend()
    else if (!is.null(colour_by)) {
        if (colour_by == "orig.ident") gg <- gg + labs(colour = "Sample")
    }

    if (!is.null(title)) gg <- gg + ggtitle(title)
    if (hide_ticks) gg <- gg + noTicks()

    return(gg)

}







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
#' @param label Logical, whether to label clusters on the plot. Default: TRUE.
#' @param palette String or character vector. If a string,
#' one of "viridis", "blues", or "redgrey", specifying which gradient
#' palette to use. Otherwise, a character vector of colours (from low to high)
#' to interpolate to create the scael. Default: redgrey.
#' @param title (Optional) String specifying the plot title
#' @param alpha Numeric, fixed alpha for points. Default: 0.6
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param legend Logical, whether or not to plot legend. Default: TRUE
#' @param hide_ticks Logical, whether to hide axis ticks. Default: FALSE
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at seurat@@ident) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#'
#' @export
#' @return A ggplot object
#'
#' @author Selin Jessa
# @examples
# tsneByMeanMarkerExpression(pbmc, "IL32")
# tsneByMeanMarkerExpression(pbmc, c("IL32", "CD2"), reduction = "pca")
tsneByMeanMarkerExpression <- function(seurat, genes,
                                       reduction = "tsne",
                                       label = TRUE,
                                       palette = "redgrey",
                                       title = NULL,
                                       alpha = 0.6,
                                       label_repel = TRUE,
                                       label_size = 4,
                                       legend = TRUE,
                                       hide_ticks = FALSE,
                                       label_short = FALSE) {

    # Get mean expression for markers
    exp_df <- meanMarkerExpression(seurat, genes)

    # Get dimensionality reduction coordinates
    exp_df <- addEmbedding(seurat, exp_df, reduction) %>%
        # Order in which points will be plot, "front" points at the bottom
        dplyr::arrange(Mean_marker_expression)

    # Get the variable names
    vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(1, 2)]

    # Plot
    gg <- exp_df %>%
        ggplot(aes(x = exp_df[[vars[1]]], y = exp_df[[vars[2]]])) +
        geom_point(aes(colour = Mean_marker_expression), size = rel(0.8), alpha = alpha)

    if (length(palette) == 1) {

        if (palette == "viridis") {

            gg <- gg + viridis::scale_color_viridis()

        } else if (palette == "blues") {

            gg <- gg + scale_colour_gradientn(
                colours = RColorBrewer::brewer.pal(n = 8, name = "Blues"))

        } else if (palette == "redgrey") {

            # NOTE: palette chosen is not the default gradient from gray -> red
            # but sets a midpoint at a lighter colour
            gg <- gg + scale_color_gradientn(
                colours = grDevices::colorRampPalette(c("gray83", "#E09797", "red"))(n = 200))

        } else {

            stop("Please pass the palette as a character vector ",
                 "or specify one of: viridis, blues, redgrey")

        }

    } else if (length(palette) == 2) {

        gg <- gg + scale_color_gradient(low = palette[1], high = palette[2])

    } else {

        gg <- gg + scale_color_gradientn(colours = palette)

    }

    if (label) {

        if (reduction == "pca") {

            message("Plotting labels is currently only available for reduction = 'tsne';",
                    " returning plot without labels.")
            return(gg)
        }

        centers <- clusterCenters(seurat)
        gg <- gg + addLabels(centers, label_repel, label_size, label_short)

    }

    axes <- gsub("_", " ", vars)

    gg <- gg +
        xlab(axes[1]) +
        ylab(axes[2]) +
        labs(colour = "Expression") +
        theme_min()

    if (!is.null(title)) gg <- gg + ggtitle(title)
    if (!legend) gg <- gg + noLegend()
    if (hide_ticks) gg <- gg + noTicks()

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
#' @param title (Optional) string used as title for the plot.
#' @param palette String or character vector. If a string,
#' one of "viridis", "blues", or "redgrey", specifying which gradient
#' palette to use. Otherwise, a character vector of colours (from low to high)
#' to interpolate to create the scael. Default: "blues".
#' @param extra Logical, plot a detailed three-panel plot, where the first
#' is a proportional bar plot of cells in each cluster in each percentile
#' group, the second is a ridge plot showing density in each cluster of the mean
#' expression of the markers (coloured by median percentile group within the
#' cluster), and the third is the labelled tSNE plot coloured
#' by percentiles. Requires \code{label = TRUE}. Default: FALSE.
#' Default: FALSE.
#' @param verbose Logical, whether to print status updates. Default: FALSE.
#' @param alpha Logical, whether to vary the alpha (point opacity) with percentile
#' group to highlight cells in the top percentiles.
#' If FALSE, sets a fixed opacity of 0.8. Default: TRUE.
#' @param legend Logical, whether to plot the legend. Default: FALSE.
#' @param point_size Numeric, size of points in scatterplot. Default: 1. (A smaller
#' value around 0.5 is better for plots which will be viewed at small scale.)
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param hide_ticks Logical, whether to hide axis ticks. Default: FALSE
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at seurat@@ident) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#'
#' @export
#'
#' @return A ggplot object
#'
#' @author Selin Jessa
#' @aliases dashboard feature
#' @examples
#' tsneByPercentileMarkerExpression(pbmc, "IL32")
#' tsneByPercentileMarkerExpression(pbmc, c("IL32", "CD2"), reduction = "pca")
#' dashboard(pbmc, "IL32", "Test dashboard")
#' feature(pbmc, "IL32")
tsneByPercentileMarkerExpression <- function(seurat, genes,
                                             label = TRUE,
                                             reduction = "tsne",
                                             palette = "blues",
                                             title = NULL,
                                             alpha = TRUE,
                                             legend = TRUE,
                                             point_size = 1,
                                             label_repel = TRUE,
                                             label_size = 4,
                                             extra = FALSE,
                                             verbose = FALSE,
                                             hide_ticks = FALSE,
                                             label_short = FALSE) {

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

    # TODO you can just set the alpha group to the gradient group,
    # and as its values, pass the discrete gradient by hand... saves some code

    # Get the variable names, and embedding
    vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(1, 2)]
    percentiles <- addEmbedding(seurat, percentiles, reduction)

    if (verbose) message("Constructing percentiles plot...")

    # Plot
    gg <- ggplot(percentiles,
                 aes(x = percentiles[[vars[1]]], y = percentiles[[vars[2]]]))


    if (alpha) {

        alpha_grad <- c(0.025, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
        alpha_grad_labels <- color_grad_labels

        gg <- gg +
            geom_point(aes(colour = factor(Gradient_group, levels = seq(1, 9)),
                           alpha = factor(Gradient_group, levels = seq(1, 9))), size = point_size) +
            scale_alpha_manual(values = alpha_grad, labels = alpha_grad_labels,
                               name = "Expression level percentile", drop = FALSE)

    } else {

        gg <- gg +
            geom_point(aes(colour = factor(Gradient_group, levels = seq(1, 9))),
                       size = point_size, alpha = 0.8) # Some fixed alpha

    }

    if (length(palette) == 1) {

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

        } else if (palette == "redgrey") {

            gg <- gg + scale_color_manual(
                values = grDevices::colorRampPalette(c("gray83", "red"))(n = 9),
                name = "Expression level percentile",
                labels = color_grad_labels,
                drop = FALSE)

        } else {

            stop("Please pass the palette as a character vector ",
                 "or specify one of: viridis, blues, redgrey")

        }

    } else {

        gg <- gg + scale_color_manual(
            values = palette,
            name = "Expression level percentile",
            labels = color_grad_labels,
            drop = FALSE)

    }

    axes <- gsub("_", " ", vars)
    gg <- gg + guides(colour = guide_legend(order = 1)) +
        xlab(axes[1]) + ylab(axes[2]) +
        theme_min()

    if (label) {

        if (reduction == "pca") {

            message("Plotting labels is currently only available for reduction = 'tsne';",
                    " returning plot without labels.")
            return(gg)
        }

        centers <- clusterCenters(seurat)
        gg <- gg + addLabels(centers, label_repel, label_size, label_short)

        if (extra) {

            # Frequency plot
            clusters <- data.frame(Cell = names(seurat@ident), Cluster = seurat@ident)
            df <- full_join(percentiles, clusters, by = "Cell")

            if (verbose) message("Computing frequencies...")

            # Calculate the proportion of cells in each cluster in each percentile group
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
                # Complete the dataframe, setting proportions to 0 when there are no cells in the cluster in a percentile group
                tidyr::complete(Cluster, Gradient_group, fill = list(n = 0, Proportion = 0, n_per_clust = 0)) %>%
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

            combined <- cowplot::plot_grid(p1, p2, gg, rel_widths = c(0.3, 0.2, 1), nrow = 1)

            if (!is.null(title)) {

                plot_title <- cowplot::ggdraw() + cowplot::draw_label(title)
                combined <- cowplot::plot_grid(plot_title, combined, ncol = 1, rel_heights = c(0.07, 1))

            }

            return(combined)

        }
    }

    if (!legend) gg <- gg + noLegend()
    if (!is.null(title)) gg <- gg + ggtitle(title)
    if (hide_ticks) gg <- gg + noTicks()

    return(gg)

}


#' @export
dashboard <- function(seurat, genes,
                      title = NULL,
                      verbose = FALSE) {


    tsneByPercentileMarkerExpression(seurat, genes, title = title,
                                     extra = TRUE, verbose = verbose)

}

#' @export
feature <- function(seurat, genes,
                    per_gene = TRUE,
                    statistic = "percentiles",
                    label = TRUE,
                    palette = "redgrey",
                    label_repel = FALSE,
                    label_size = 4,
                    label_short = FALSE,
                    legend = FALSE,
                    title = NULL,
                    reduction = "tsne",
                    alpha = ifelse(statistic == "percentiles", FALSE, 0.6),
                    point_size = 0.5,
                    ncol = ifelse(length(genes) == 1, 1, ifelse(length(genes) %in% c(2, 4), 2, 3)),
                    hide_ticks = FALSE) {

    if ((length(genes) >= 20) & per_gene) message("NOTE: you have input a lot of genes! ",
                                                  "This function by default generates ",
                                                  "one plot per gene. Set per_gene = FALSE ",
                                                  "to plot a summary statistic of all genes.")

    if (statistic == "percentiles")  {

        if (per_gene) {

            genes_out <- findGenes(seurat, genes)
            if (length(genes_out$undetected > 0)) message(paste0("NOTE: [",
                                                                 paste0(genes_out$undetected, collapse = ", "),
                                                                 "] undetected in the data"))

            if(length(genes_out$detected) == 0) stop("No genes specified were ",
                                                     "found in the data.")

            if ((ncol == 3) & (length(genes_out$detected) < 3)) ncol <- 2

            plots <- cowplot::plot_grid(
                plotlist = lapply(genes_out$detected,
                                  function(gene) tsneByPercentileMarkerExpression(seurat,
                                                                                  gene,
                                                                                  label = label,
                                                                                  palette = palette,
                                                                                  label_repel = label_repel,
                                                                                  label_size = label_size,
                                                                                  label_short = label_short,
                                                                                  title = gene,
                                                                                  alpha = alpha,
                                                                                  legend = legend,
                                                                                  point_size = point_size,
                                                                                  hide_ticks = hide_ticks)),
                ncol = ncol)

            if (is.null(title)) return(plots)
            else {

                plot_title <- cowplot::ggdraw() + cowplot::draw_label(title, hjust = 0, size = 12)
                cowplot::plot_grid(plot_title, plots, ncol = 1, rel_heights = c(0.05, 1))

            }


        } else {

            tsneByPercentileMarkerExpression(seurat,
                                             genes,
                                             label = label,
                                             palette = palette,
                                             label_repel = label_repel,
                                             label_size = label_size,
                                             label_short = label_short,
                                             title = title,
                                             alpha = alpha,
                                             legend = legend,
                                             point_size = point_size,
                                             hide_ticks = hide_ticks)

        }


    } else if (statistic == "mean") {

        if (per_gene) {

            genes_out <- findGenes(seurat, genes)
            if (length(genes_out$undetected > 0)) message(paste0("NOTE: [",
                                                                 paste0(genes_out$undetected, collapse = ", "),
                                                                 "] undetected in the data"))

            if(length(genes_out$detected) == 0) stop("No genes specified were ",
                                                     "found in the data.")

            if ((ncol == 3) & (length(genes_out$detected) < 3)) ncol <- 2

            plots <- cowplot::plot_grid(
                plotlist = lapply(genes_out$detected,
                                  function(gene) tsneByMeanMarkerExpression(seurat,
                                                                            gene,
                                                                            reduction = reduction,
                                                                            palette = palette,
                                                                            title = gene,
                                                                            legend = legend,
                                                                            label = label,
                                                                            label_short = label_short,
                                                                            label_repel = label_repel,
                                                                            label_size = label_size,
                                                                            hide_ticks = hide_ticks)),
                ncol = ncol)

            if (is.null(title)) return(plots)
            else {

                plot_title <- cowplot::ggdraw() + cowplot::draw_label(title, hjust = 0, size = 12)
                cowplot::plot_grid(plot_title, plots, ncol = 1, rel_heights = c(0.05, 1))

            }

        } else {

            tsneByMeanMarkerExpression(seurat,
                                       genes,
                                       reduction = reduction,
                                       palette = palette,
                                       title = title,
                                       legend = legend,
                                       label = label,
                                       label_repel = label_repel,
                                       label_size = label_size,
                                       label_short = label_short,
                                       hide_ticks = hide_ticks)
        }
    }
}


#' vln
#'
#' Similar to Seurat::VlnPlot() except it prints genes that are not found in the
#' data, but continues plotting without error,
#' and has a (default) option to group plots by cluster instead of gene.
#'
#' @param seurat Seurat object
#' @param genes Genes to plot violins for
#' @param facet_by String, one of "gene" or "cluster". Default: "cluster",
#' genes will be on the x axis, and there will be one plot per cluster.
#' If "gene", clusters will be on the x axis, and there will be one plot per
#' gene (akin to Seurat::VlnPlot)
#' @param point_size Numeric value for point size, use -1 to hide points. Default: 0.1.
#' @param adjust Bandwidth for density estimation, passed to \code{geom_violin}.
#' See ggplot2 documentation for more info. Default: 1. NOTE/TODO: If vln() with
#' facet = gene looks different from Seurat::VlnPlot, this is probably the culprit.
#'
#' @return A ggplot2 object
#' @author Selin Jessa
#' @export
#'
#' @examples
#' vln(pbmc, c("IL32", "MS4A1"))
#' vln(pbmc, c("IL32", "MS4A1"), facet_by = "gene")
#' vln(pbmc, c("IL32", "MS4A1"), point_size = -1, facet_by = "gene")
vln <- function(seurat, genes, facet_by = "cluster", point_size = 0.1, adjust = 1) {

    exp <- fetchData(seurat, genes, return_cluster = TRUE) %>%
        tidyr::gather(Gene, Expression, 2:length(.))

    genes_out <- findGenes(seurat, genes)

    if (facet_by == "gene") gg <- ggplot(exp, aes(x = Cluster, y = Expression, fill = Cluster))
    else if (facet_by == "cluster") gg <- ggplot(exp, aes(x = Gene, y = Expression, fill = Cluster))

    # Seurat::SingleVlnPlot limits the data in this way
    # y.max <- max(select(exp, Expression))
    # y.min <- min(select(exp, Expression))

    gg <- gg +
        geom_violin(scale = "width", adjust = adjust, trim = TRUE) +
        geom_jitter(size = point_size, alpha = 0.5) +
        theme_min() +
        theme(legend.position = "none")
    # + ylim(y.min, y.max)

    if (facet_by == "gene") {

        gg <- gg + facet_wrap(~ Gene, ncol = 2) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

    } else if (facet_by == "cluster") gg <- gg + facet_wrap(~ Cluster, scales = "free_y")

    gg <- gg + rotateX()
    return(gg)

}



#' vlnGrid
#'
#' A grid of small violin plots, one plot per gene per cluster, similar to
#' Figure 5D in the Drop-seq paper by Macosko et al.
#'
#' @param seurat Seurat object
#' @param genes Genes to plot violins for, in the order in which they should be
#' plotted.
#' @param order Either "genes", or a vector containing the clusters in the order
#' they should be plotted (from top to bottom). If "genes", clusters will be sorted
#' in decreasing order of their median expression* of these genes, in the order
#' in which they're provided. Passing a vector will cause clusters to be plot
#' in that order. Default: "genes". *We take the median of cells with non-zero
#' expression values.
#' @param subset_clusters (Optional) Vector of clusters to include on the plot.
#' @param width String, one of "width", "area", or "count". From the ggplot2
#' documentation of \code{geom_violin}: "if "area" (default), all violins have
#' the same area (before trimming the tails). If "count", areas are scaled
#' proportionally to the number of observations. If "width", all violins have
#' the same maximum width." Default: "width".
#'
#' @return A ggplot2 object
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # Use the first 15 genes in the data, order by gene
#' vlnGrid(pbmc, head(rownames(pbmc@data), 15))
#'
#' # Specify order
#' vlnGrid(pbmc, head(rownames(pbmc@data), 15), order = c(0, 3, 2, 1))
#'
#' # Plot only a subset of clusters
#' vlnGrid(pbmc, head(rownames(pbmc@data), 15), subset_clusters = c(0, 1))
vlnGrid <- function(seurat, genes,
                    subset_clusters = NULL,
                    order = "genes",
                    scale = "width",
                    title = NULL) {

    expr <- cytokit::fetchData(seurat, genes, return_cell = TRUE, return_cluster = TRUE) %>%
        tidyr::gather(Marker, Expression, 3:length(.))

    # return(expr)

    if (!is.null(subset_clusters)) expr <- filter(expr, Cluster %in% subset_clusters)

    genes_out <- findGenes(seurat, genes)

    if (length(order) == 1) {

        if (order == "genes") {

            sort_criteria <- glue::glue("desc({genes_out$detected})")

            cluster_order <- expr %>%
                filter(Expression != 0) %>%
                tidyr::spread(Marker, Expression) %>%
                select(-Cell) %>%
                group_by(Cluster) %>%
                summarise_all(funs(median), na.rm = TRUE) %>%
                arrange_(sort_criteria) %>%
                `[[`("Cluster")

            expr$Cluster <- factor(expr$Cluster, levels = rev(cluster_order))

        } else stop("Please set order = 'genes', or provide an ordering of clusters.")

    } else {

        if ((is.null(subset_clusters)) && (length(order) != length(unique(seurat@ident)))) {

            stop("Please provide an ordering which includes all clusters.")

        } else if ((!is.null(subset_clusters)) && (length(order) != length(subset_clusters))) {

            stop("Please provide an ordering which includes all clusters passed to ",
                 "the 'subset_clusters' argument.")
        }

        expr$Cluster <- factor(expr$Cluster, levels = rev(order))

    }

    expr$Marker <- factor(expr$Marker, levels = genes_out$detected)

    # Set colours to the levels of the clusters, so that they are preserved
    # even if an order was specified
    colours <- ggColours(length(levels(seurat@ident)))
    names(colours) <- levels(seurat@ident)

    gg <- expr %>%
        ggplot(aes(x = Cluster, y = Expression)) +
        geom_violin(aes(fill = Cluster), scale = scale, size = 0.5) +
        scale_fill_manual(values = colours) +
        facet_wrap(~ Marker, ncol = length(unique(expr$Marker))) +
        theme_min() +
        coord_flip() +
        ggplot2::theme(panel.border = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.text.x = element_blank(),
                       axis.line.x = element_blank(),
                       legend.position = "none",
                       strip.text.x = element_text(angle = 30, size = 8))

    if (!is.null(title)) gg <- gg + ggtitle(title)

    return(gg)

}





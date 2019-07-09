# Functions for plotting the data with computed information



#' Colour cells in t-SNE or PCA space by gene expression
#'
#' Plot a low-dimensional embedding of the cells,
#' coloured by expression of a gene, or mean expression of a group of marker
#' genes. Defaults to t-SNE space, but see the \code{reduction} argument for
#' how to plot in PCA space instead. This function is based on \code{Seurat::FeaturePlot}.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunPCA() or Seurat::RunTSNE() to the object)
#' @param genes String or character vector specifying gene(s) to use
#' @param reduction String specifying the dimensionality reduction to use,
#' retrieves t-SNE by default. This should match the names of the elements of
#' the list seurat@@dr, so it will typically be one of "pca" or "tsne".
#' Default: "tsne"
#' @param dim1 Numeric, index of dimension from \code{reduction} to plot on
#' the x-axis. e.g. to plot the 3rd PC on the x-axis, pass 3. Default: 1.
#' @param dim2 Numeric, like \code{dim2}, but for the y-axis. Default: 2.
#' @param label Logical, whether to label clusters on the plot. Default: TRUE.
#' @param palette String or character vector. If a string,
#' one of "viridis", "blues", or "redgrey", specifying which gradient
#' palette to use. Otherwise, a character vector of colours (from low to high)
#' to interpolate to create the scael. Default: redgrey.
#' @param title (Optional) String specifying the plot title
#' @param alpha Numeric, fixed alpha for points. Default: 0.6
#' @param point_size Numeric, size of points in scatterplot. Default: 1. (A smaller
#' value around 0.5 is better for plots which will be viewed at small scale.)
#' @param label_repel Logical, if \code{label} is TRUE, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param legend Logical, whether or not to plot legend. Default: TRUE
#' @param hide_ticks Logical, whether to hide axis ticks. Default: FALSE
#' @param hide_axes Logical, whether to hide axis labels. Default: TRUE
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at seurat@@ident) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param limits (Optional) A numeric vector of length two providing the limits to
#' use for the colour scale (documentation
#' from \code{\link[ggplot2]{continous_scale}}. Default: 0 and max of the data.
#'
#' @export
#' @return A ggplot object
#'
#' @author Selin Jessa
#' @examples
#' tsneByMeanMarkerExpression(pbmc, "IL32")
#' tsneByMeanMarkerExpression(pbmc, c("IL32", "CD2"), reduction = "pca")
#' tsneByMeanMarkerExpression(pbmc, "IL32", reduction = "pca", dim1 = 1, dim2 = 3)
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
                                       hide_axes = FALSE,
                                       limits = NULL,
                                       label_short = FALSE,
                                       dim1 = 1,
                                       dim2 = 2,
                                       return_df = FALSE,
                                       point_size = 1) {

    # Get mean expression for markers
    exp_df <- meanGeneExpression(seurat, genes)

    # Get dimensionality reduction coordinates
    exp_df <- addEmbedding(seurat, exp_df, reduction, dim1, dim2) %>%
        # Order in which points will be plot, "front" points at the bottom
        dplyr::arrange(Mean_marker_expression)

    if (return_df) return(exp_df)

    # Get the variable names
    vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(dim1, dim2)]

    # Set limits: if not provided, use default min/max
    if (is.null(limits)) limits <- c(NA, NA)

    # Plot
    gg <- exp_df %>%
        ggplot(aes(x = exp_df[[vars[1]]], y = exp_df[[vars[2]]])) +
        geom_point(aes(colour = Mean_marker_expression), size = point_size, alpha = alpha)

    if (length(palette) == 1) {

        if (palette == "viridis") {

            gg <- gg + viridis::scale_color_viridis(limits = limits)

        } else if (palette == "blues") {

            gg <- gg + scale_colour_gradientn(
                colours = RColorBrewer::brewer.pal(n = 8, name = "Blues"),
                limits = limits)

        } else if (palette == "redgrey") {

            # NOTE: palette chosen is not the default gradient from gray -> red
            # but sets a midpoint at a lighter colour
            gg <- gg + scale_color_gradientn(
                colours = grDevices::colorRampPalette(c("gray83", "#E09797", "red"))(n = 200),
                limits = limits)

        } else {

            stop("Please pass the palette as a character vector ",
                 "or specify one of: viridis, blues, redgrey")

        }

    } else if (length(palette) == 2) {

        gg <- gg + scale_color_gradient(low = palette[1], high = palette[2], limits = limits)

    } else {

        gg <- gg + scale_color_gradientn(colours = palette, limits = limits)

    }


    if (label) {

        centers <- clusterCenters(seurat, reduction = reduction, dim1 = dim1, dim2 = dim2)
            gg <- gg + addLabels(centers, label_repel, label_size, label_short)

    }

    if (reduction == "tsne") axes <- gsub("_", " ", vars)
    else if (reduction == "pca") {

        pve <- getVarianceExplained(seurat) %>%
            .$percent.var.explained %>%
            round(1)

        axes <- c(glue('{gsub("_", " ", vars[1])} ({pve[dim1]}%)'),
                  glue('{gsub("_", " ", vars[2])} ({pve[dim2]}%)'))

    }

    gg <- gg +
        xlab(axes[1]) +
        ylab(axes[2]) +
        labs(colour = "Expression") +
        theme_min()

    if (!is.null(title)) gg <- gg + ggtitle(title)
    if (!legend) gg <- gg + noLegend()
    if (hide_ticks) gg <- gg + noTicks()
    if (hide_axes) gg <- gg + xlab(NULL) + ylab(NULL)

    return(gg)

}


#' Colour cells in t-SNE space by percentile gene expression
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
#' @param legend_options String, "percentiles" or "values". Default: "percentiles".
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
#' @aliases dashboard
#' @examples
#' tsneByPercentileMarkerExpression(pbmc, "IL32")
#' dashboard(pbmc, "IL32", title = "Test dashboard")
#' feature(pbmc, "IL32")
tsneByPercentileMarkerExpression <- function(seurat, genes,
                                             reduction = "tsne",
                                             label = TRUE,
                                             palette = "blues",
                                             title = NULL,
                                             alpha = TRUE,
                                             legend = TRUE,
                                             legend_options = "percentiles",
                                             point_size = 1,
                                             label_repel = TRUE,
                                             label_size = 4,
                                             extra = FALSE,
                                             verbose = FALSE,
                                             hide_ticks = FALSE,
                                             label_short = FALSE,
                                             dim1 = 1,
                                             dim2 = 2) {

    if (verbose) message("Computing percentiles...")

    # Get expression percentiles
    percentiles <- percentilesMarkerExpression(seurat, genes)

    if (legend_options=="percentiles") {
        color_grad_labels <- c("Undetected",
                               "> 0 & \\u2264 50",
                               "> 50 & \\u2264 70",
                               "> 70 & \\u2264 90",
                               "> 90 & \\u2264 92",
                               "> 92 & \\u2264 94",
                               "> 94 & \\u2264 96",
                               "> 96 & \\u2264 98",
                               "> 98 & \\u2264 100")
    } else if (legend_options=="values") {
	# If legend_options is set to values, compute the values corresponding to the percentiles.
	# Cell.type corresponds to the values. This variable should probably be renamed.
        labels.min <- group_by(percentiles, Gradient_group) %>% summarize(minValue=min(Cell.type)) %>% .$minValue
        labels.max <- group_by(percentiles, Gradient_group) %>% summarize(maxValue=max(Cell.type)) %>% .$maxValue
        color_grad_labels <- c("Undetected", paste0("> ", round(labels.min[2:length(labels.min)], digits=2), " & \u2264 ", round(labels.max[2:length(labels.max)], digits=2)))
    }

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
                               name = "Expression level percentile",
                               drop = ifelse(legend_options == "percentiles", FALSE, TRUE))

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
                    drop = ifelse(legend_options == "percentiles", FALSE, TRUE))

        } else if (palette == "blues") {

            gg <- gg +
                scale_color_manual(
                    values = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
                    name = "Expression level percentile",
                    labels = color_grad_labels,
                    drop = ifelse(legend_options == "percentiles", FALSE, TRUE))

        } else if (palette == "redgrey") {

            gg <- gg + scale_color_manual(
                values = grDevices::colorRampPalette(c("gray83", "red"))(n = 9),
                name = "Expression level percentile",
                labels = color_grad_labels,
                drop = ifelse(legend_options == "percentiles", FALSE, TRUE))

        } else {

            stop("Please pass the palette as a character vector ",
                 "or specify one of: viridis, blues, redgrey")

        }

    } else {

        gg <- gg + scale_color_manual(
            values = palette,
            name = "Expression level percentile",
            labels = color_grad_labels,
            drop = ifelse(legend_options == "percentiles", FALSE, TRUE))

    }

    axes <- gsub("_", " ", vars)
    gg <- gg + guides(colour = guide_legend(order = 1)) +
        xlab(axes[1]) + ylab(axes[2]) +
        theme_min()

    if (label) {

        centers <- clusterCenters(seurat, reduction = reduction, dim1 = dim1, dim2 = dim2)
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
            df_means <- meanGeneExpression(seurat, genes) %>%
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

            combined <- plot_grid(p1, p2, gg, rel_widths = c(0.3, 0.2, 1), nrow = 1)

            if (!is.null(title)) {

                plot_title <- cowplot::ggdraw() + cowplot::draw_label(title)
                combined <- plot_grid(plot_title, combined, ncol = 1, rel_heights = c(0.07, 1))

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
dashboard <- function(seurat,
                      genes,
                      palette = "viridis",
                      title = NULL,
                      verbose = FALSE) {


    tsneByPercentileMarkerExpression(seurat, genes,
                                     palette = palette,
                                     title = title,
                                     extra = TRUE, verbose = verbose)

}

#' @describeIn tsneByMeanMarkerExpression Shortcut function for plotting mean expression
#' @export
feature <- function(seurat, genes,
                    per_gene = TRUE,
                    statistic = "mean",
                    label = TRUE,
                    palette = "redgrey",
                    label_repel = FALSE,
                    label_size = 4,
                    label_short = FALSE,
                    legend = FALSE,
                    title = NULL,
                    reduction = "tsne",
                    limits = c(NA, NA),
                    dim1 = 1,
                    dim2 = 2,
                    alpha = ifelse(statistic == "percentiles", FALSE, 0.6),
                    point_size = 0.5,
                    ncol = ifelse(length(genes) == 1, 1, ifelse(length(genes) %in% c(2, 4), 2, 3)),
                    hide_ticks = TRUE,
                    hide_axes = TRUE) {

    if ((length(genes) >= 20) & per_gene) message("NOTE: you have input a lot of genes! ",
                                                  "This function by default generates ",
                                                  "one plot per gene. Set per_gene = FALSE ",
                                                  "to plot a summary statistic of all genes.")

    if (statistic == "percentiles")  {

        if (reduction != "pca") warning("Mapping expression by percentiles ",
                                        "is not yet implemented for ")

        if (per_gene) {

            genes_out <- findGenes(seurat, genes)
            if (length(genes_out$undetected > 0)) message(paste0("NOTE: [",
                                                                 paste0(genes_out$undetected, collapse = ", "),
                                                                 "] undetected in the data"))

            if(length(genes_out$detected) == 0) stop("No genes specified were ",
                                                     "found in the data.")

            if ((ncol == 3) & (length(genes_out$detected) < 3)) ncol <- 2

            plots <- plot_grid(
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
                plot_grid(plot_title, plots, ncol = 1, rel_heights = c(0.05, 1))

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

            plots <- plot_grid(
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
                                                                            hide_ticks = hide_ticks,
                                                                            hide_axes = hide_axes,
                                                                            point_size = point_size,
                                                                            limits = limits,
                                                                            dim1 = dim1,
                                                                            dim2 = dim2)),
                ncol = ncol)

            if (is.null(title)) return(plots)
            else {

                plot_title <- cowplot::ggdraw() + cowplot::draw_label(title, hjust = 0, size = 12)
                plot_grid(plot_title, plots, ncol = 1, rel_heights = c(0.05, 1))

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
                                       hide_ticks = hide_ticks,
                                       hide_axes = hide_axes,
                                       point_size = point_size,
                                       limits = limits,
                                       dim1 = dim1,
                                       dim2 = dim2)
        }
    }
}


#' Generate violin plots of gene expression in each cluster
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



#' Generate a grid of tiny violins of gene expression in each cluster
#'
#' A grid of small violin plots, one plot per gene per cluster, similar to
#' Figure 5D in the Drop-seq paper by Macosko et al. Aka "fish plot".
#'
#' @param seurat Seurat object
#' @param genes Genes to plot violins for, in the order (left to right)
#' in which they should be  plotted.
#' @param order Either "genes", or a vector containing the clusters in the order
#' they should be plotted (from top to bottom). If "genes", clusters will be sorted
#' in decreasing order of their median expression* of these genes, in the order
#' in which they're provided. Passing a vector will cause clusters to be plot
#' in that order. Default: "genes". *We take the median of cells with non-zero
#' expression values (in \code{seurat@@data}).
#' @param subset_clusters (Optional) Vector of clusters to include on the plot.
#' @param colours (Optional) Vector of colours to use. Either one colour
#' per cluster, in the order of \code{levels(seurat@@ident)}, or one colour per
#' cluster passed to \code{subset_clusters}, in the other they were provided.
#' Default: default ggplot2 colours used by Seurat.
#' @param width String, one of "width", "area", or "count". From the ggplot2
#' documentation of \code{\link[ggplot2]{geom_violin}}: "if "area" (default), all violins have
#' the same area (before trimming the tails). If "count", areas are scaled
#' proportionally to the number of observations. If "width", all violins have
#' the same maximum width." Default: "width".
#' @param scales See \code{scales} parameter in \code{\link[ggplot2]{facet_wrap}}.
#' Default: scales = "free_x", which allows every gene to have its own scale.
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
#'
#' # Specify colours
#' vlnGrid(pbmc, head(rownames(pbmc@data), 15), subset_clusters = c(0, 1),
#'         colours = c("red", "blue"))
vlnGrid <- function(seurat, genes,
                    subset_clusters = NULL,
                    order = "genes",
                    colours = NULL,
                    scale = "width",
                    title = NULL,
                    scales = "free_x") {

    expr <- cytobox::fetchData(seurat, genes, return_cell = TRUE, return_cluster = TRUE) %>%
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


    if (!is.null(colours) && length(colours) != length(levels(seurat@ident))) {

        if (length(colours) != length(subset_clusters)) {

            stop("Not enough colours! Please provide as many colours ",
                 "as clusters in the dataset, or one per cluster specified in the ",
                 "'subset_clusters' argument.")

        }

    } else if (is.null(colours)) {

        colours <- ggColours(length(levels(seurat@ident)))
        names(colours) <- levels(seurat@ident)

    }

    gg <- expr %>%
        ggplot(aes(x = Cluster, y = Expression)) +
        geom_violin(aes(fill = Cluster), scale = scale, size = 0.5) +
        scale_fill_manual(values = colours) +
        facet_wrap(~ Marker, ncol = length(unique(expr$Marker)),
                   scales = scales) +
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





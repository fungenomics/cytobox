# Functions for basic visualization of the data



#' Plot a tSNE embedding for a dataset
#'
#' Plot a dataset in tSNE space, akin to Seurat::TSNEPlot()
#'
#' @param seurat Seurat object, where Seurat::RunTSNE() has been applied
#' @param colour_by (Optional) String, specifying the column in \code{seurat@@meta.data}
#' by which to colour cells. Default: NULL, colour cells by cluster (in \code{seurat@@ident}).
#' @param colours (Optional) Character vector of colours for points. If \code{colour_by}
#' is NULL, cells will be coloured by cluster; this should be a named character vector of colours for points. Names should
#' correspond to cluster names (e.g. \code{levels(seurat@@ident)}). If
#' specifying \code{colour_by}, and the variable is discrete, this should be a character vector,
#' with either names or order corresponding to categorical
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
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes cluster
#' names (at seurat@@ident) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param cells Character vector of cell names if only a subset of cells should be
#' plot (these should correspond to seurat@@cell.names). Default: Plot all cells.
#' See the argument \code{clusters_to_label} for only labelling certain clusters.
#' See the \code{constrain_scale} argument for controlling the scales of the plot.
#' @param order_by String, corresponding to a column in seurat@@meta.data, specifying
#' a variable to control the order in which cells are plot. (Thus, you can manually
#' specify the order, add it as a new column in seurat@@meta.data, and pass that).
#' If numeric, cells with high values are plot on top. If not, the column must
#' be a factor, and cells will be ordered according to the levels, with cells
#' in the first level plot on top. Default: if a numeric column is specified
#' to \code{colour_by}, sort by that variable, otherwise, use the ordering of the cells
#' in the Seurat object.
#' @param clusters_to_label (Optional.) If \code{label} is TRUE,
#' clusters for which labels should be plot (if only a subset of clusters should be labelled).
#' Default: NULL (Label all clusters).
#' @param na_color String, specifying the colour (built-in or hex code) to use to
#' plot points which have an NA value, for example
#' in the variable specified in \code{colour_by}. Default: light gray ("gray80),
#' change to "white" to purposely hide those cells. If you do not want to plot
#' certain cells at all, pass names of cells to plot to the \code{cells} argument.
#' @param limits Numeric vector of length two providing the lower and upper limits of
#' the colour scale, if colouring by a continuous variable. Default: min and max
#' of the values the variable takes on in the data.
#' @param constrain_scale Logical, if plotting a subset of cells, whether to
#' use the limits of the tSNE embedding computed on the whole dataset (useful
#' for constraining scales across plots while only plotting specific cells).
#' Default: TRUE
#'
#' @return A ggplot2 object
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' tsne(pbmc)
#'
#' # Demonstrate label_short:
#' # Set cluster IDs to be too long
#' pbmc2 <- pbmc
#' levels(pbmc2@ident) <- c("1-Cell type A", "2-Cell type B", "3-Cell type C", "4-Cell type D")
#' tsne(pbmc2)
#'
#' # Plot the prefixes only
#' tsne(pbmc2, label_short = TRUE)
tsne <- function(seurat,
                 colour_by = NULL,
                 colours = NULL,
                 colour_by_type = "discrete",
                 label = TRUE,
                 point_size = 0.6, alpha = 0.8,
                 legend = ifelse((is.null(colour_by)) && (label), FALSE, TRUE),
                 label_repel = TRUE,
                 label_size = 4,
                 cells = NULL,
                 order_by = NULL,
                 clusters_to_label = NULL,
                 hide_ticks = FALSE,
                 title = NULL,
                 label_short = FALSE,
                 na_color = "gray80",
                 limits = NULL,
                 constrain_scale = TRUE) {

    embedding <- data.frame(Cell = seurat@cell.names,
                            tSNE_1 = seurat@dr$tsne@cell.embeddings[, 1],
                            tSNE_2 = seurat@dr$tsne@cell.embeddings[, 2],
                            Cluster = seurat@ident,
                            stringsAsFactors = FALSE)

    if (!is.null(order_by)) {

        # Check the variable is usable for sorting
        if (!is.numeric(seurat@meta.data[[order_by]]) && !is.factor(seurat@meta.data[[order_by]])) {

            stop("The variable specified in 'order_by' is neither numeric ",
                 "nor a factor. If the column is of type character, consider ",
                 "converting it to a factor. Otherwise, pass the name of a numeric column.")

        }

        embedding[[order_by]] <- seurat@meta.data[[order_by]]

    } else if ((!is.null(colour_by)) && is.numeric(seurat@meta.data[[colour_by]])) {

        # If order_by is not specified but colour_by is, and is numeric,
        # by default, order the cells by that variable
        order_by <- colour_by

    }

    if (!is.null(colour_by)) embedding[[colour_by]] <- seurat@meta.data[[colour_by]]
    if (!is.null(cells)) embedding <- embedding %>% filter(Cell %in% cells)
    if (!is.null(order_by)) embedding <- embedding %>% arrange_(order_by)

    gg <- ggplot(embedding, aes(x = tSNE_1, y = tSNE_2))

    if (label && all(is.na(seurat@ident))) {
        label <- FALSE
        message("NOTE: identity of all cells is NA, setting 'label' to FALSE.")
    }

    if (is.null(colour_by)) {

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

        if (is.null(limits)) lims <- c(NA, NA)
        else lims <- limits

        gg <- gg +
            geom_point(aes_string(colour = colour_by), size = point_size, alpha = alpha)

        if (!is.null(colours)) {

            if (colour_by_type == "discrete") gg <- gg + scale_color_manual(values = colours, na.value = na_color)
            else if (colour_by_type == "continuous") {

                gg <- gg + scale_color_gradientn(colours = colours,
                                                 na.value = na_color,
                                                 limits = lims)
            }

        } else {

            if (colour_by_type == "continuous") { # Otherwise for discrete, default ggplot2 colours are used

                gg <- gg + scale_color_gradientn(colours = grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "OrRd"))(n = 100),
                                                 na.value = na_color,
                                                 limits = lims)

            }
        }

    }

    if (label) {

        centers <- clusterCenters(seurat)
        gg <- gg + addLabels(centers     = centers,
                             label_repel = label_repel,
                             label_size  = label_size,
                             label_short = label_short,
                             clusters    = clusters_to_label)

    }

    gg <- gg + theme_min() + xlab("tSNE 1") + ylab("tSNE 2")

    if (!legend) gg <- gg + noLegend()
    else if (!is.null(colour_by)) {
        if (colour_by == "orig.ident") gg <- gg + labs(colour = "Sample")
    }

    if (!is.null(title)) gg <- gg + ggtitle(title)
    if (hide_ticks) gg <- gg + noTicks()
    if (constrain_scale) gg <- gg + constrainScale(seurat, reduction = "tsne")

    return(gg)

}




#' Highlight select clusters on a tSNE plot
#'
#' Plot the tSNE embedding for the dataset, with cells from select clusters coloured
#' by either their original colours or provided colours, and cells from all
#' other clusters in another (non-intrusive) colour, or not at all. This is a
#' thin wrapper for \code{\link{tsne}} which takes care of specifying cells and colours
#' in order to highlight the desired clusters.
#'
#' @param seurat Seurat object, where Seurat::RunTSNE() has been applied
#' @param clusters Vector of one or more clusters to highlight, matching the levels at
#' \code{levels(seurat@@ident)}. If "none", all clusters are coloured by \code{default_colour}.
#' @param original_colours (Optional) Vector of colours to use. Either one colour
#' per cluster, in the order of \code{levels(seurat@@ident)}, or one colour per
#' cluster passed to \code{clusters}, in the other they were provided.
#' Default: default ggplot2 colours used by Seurat.
#' @param default_colour String, colour to use for non-highlighted clusters, or
#' "none", if cells in those clusters should not be plot at all. Default: gray80
#' (light grey).
#' @param label_all Logical, if labelling the tSNE (if \code{label == TRUE}), whether
#' to label all the clusters, or only the ones being highlighted. Default: FALSE.
#'
#' @inheritDotParams tsne -seurat -clusters_to_label -colours -label
#' @return A ggplot2 object
#' @export
#' @author Selin Jessa
#'
#' @examples
#' # Highlight cluster 3 on the tSNE plot
#' highlight(pbmc, 3)
#'
#' # Pass additional arguments to cytokit::tsne
#' highlight(pbmc, c(2, 3), label = FALSE, title = "Test highlight")
#'
#' # Change default colour
#' highlight(pbmc, c(2, 3), default_colour = "lightblue")
#'
#' # Specify the colours to highlight the clusters with
#' highlight(pbmc, c(2, 3), c("red", "blue"))
#'
#' # Don't plot cells in other clusters
#' highlight(pbmc, c(2, 3), default_colour = "none")
highlight <- function(seurat,
                      clusters,
                      original_colours = NULL,
                      default_colour = "gray90",
                      label_all = FALSE,
                      label = TRUE,
                      ...) {

    highlight_colours <- getClusterColours(seurat = seurat,
                                           clusters = clusters,
                                           original_colours = original_colours,
                                           default_colour = default_colour)

    cells_label <- whichCells(seurat, clusters)

    if (label_all) {

        # Assign an order to the cells so the highlighted ones will be plot on top
        seurat@meta.data$Order_cells <- ifelse(seurat@cell.names %in% cells_label,
                                               as.character(seurat@ident),
                                               "Other")

        seurat@meta.data$Order_cells <- factor(seurat@meta.data$Order_cells, levels = c(clusters, "Other"))

        tsne(seurat, colours = highlight_colours, order_by = "Order_cells", ...)

    }

    else if (default_colour == "none") tsne(seurat,
                                            colours = highlight_colours,
                                            clusters_to_label = clusters,
                                            # Don't plot cells from the non-highlighted clusterss
                                            cells = cells_label,
                                            label = label,
                                            ...)
    else tsne(seurat,
              colours = highlight_colours,
              clusters_to_label = clusters,
              label = label,
              ...)

}



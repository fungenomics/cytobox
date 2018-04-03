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




#' Highlight select clusters on a tSNE plot
#'
#' Plot the tSNE embedding for the dataset, with cells from select clusters coloured
#' by either their original colours or provided colours, and cells from all
#' other clusters in another (non-intrusive) colour.
#'
#' @param seurat Seurat object, where Seurat::RunTSNE() has been applied
#' @param clusters Vector of one or more clusters to highlight, matching the levels at
#' \code{levels(seurat@@ident)}. If "none", all clusters are coloured by \code{default_colour}.
#' @param original_colours (Optional) Vector of colours to use. Either one colour
#' per cluster, in the order of \code{levels(seurat@@ident)}, or one colour per
#' cluster passed to \code{clusters}, in the other they were provided.
#' Default: default ggplot2 colours used by Seurat.
#' @param default_colour Colour to use for non-highlighted clusters. Default: gray80
#' (light grey).
#' @param ... (Optional) Other arguments passed to \code{\link{tsne}}.
#'
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
highlight <- function(seurat, clusters, original_colours = NULL, default_colour = "gray80", ...) {

    highlight_colours <- getClusterColours(seurat = seurat,
                                           clusters = clusters,
                                           original_colours = original_colours,
                                           default_colour = default_colour)

    cytokit::tsne(seurat, colours = highlight_colours, ...)

}



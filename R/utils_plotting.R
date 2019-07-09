# Helper and utility functions for plotting



#' Get default ggplot2/Seurat colours
#'
#' Get evenly spaced colours from around the colour wheel, which are the default
#' colours assigned to clusters by Seurat. The output of this function can be
#' passed to the \code{scale_colour_manual()} and \code{scale_fill_manual()} functions
#' from ggplot2, as the \code{values} argument. (\code{\link{ggColors}} points
#' to this function.)
#'
#' @param n Number of colours to return
#'
#' @return Named character vector, where names are the names of clusters, from
#' 0 to n-1, and values are the hex codes for the colours.
#' @export
#'
#' @examples
#'
#' n_clust <- 5
#' ggColours(n_clust)
#'
#' @references https://stackoverflow.com/a/8197703
#' @aliases ggColors
#' @importFrom grDevices hcl
ggColours <- function(n) {

    hues <- seq(15, 375, length = n + 1)
    colours <- hcl(h = hues, l = 65, c = 100)[1:n]
    names(colours) <- seq(0, n - 1) # Since the first cluster in Seurat is 0

    return(colours)

}

#' @export
ggColors <- ggColours




#' Get a vector of cluster colours, optionally highlighting select clusters
#'
#' Given a seurat object, get a named character vector of cluster colours, where
#' names are cluster names (coresponding to \code{levels(seurat@@ident)}), and
#' values are hex codes of the colours, either the default colours from Seurat,
#' or colours you specify. This is trivial to make this yourself,
#' so this is used as a utility function for retaining the colours of
#' clusters you want to highlight, and setting the colours of all other clusters
#' to grey, or another default, non-intrusive colour.
#'
#' @param seurat Seurat object
#' @param clusters Vector of one or more clusters to highlight, matching the levels at
#' \code{levels(seurat@@ident)}. If "none", returns \code{default_colour}
#' for every clusters. Default: all clusters, obtained from \code{levels(seurat@@ident)}.
#' @param original_colours  (Optional) Vector of colours to use. Either one colour
#' per cluster, in the order of \code{levels(seurat@@ident)}, or one colour per
#' cluster passed to \code{clusters}, in the other they were provided.
#' Default: default ggplot2 colours used by Seurat.
#' @param default_colour Colour to use for non-highlighted clusters
#' Default: gray80 (light grey).
#'
#' @return Named character vector
#' @export
#' @author Selin Jessa
#' @examples
#'
#' # Trivial: get named character vector with default colours
#' getClusterColours(pbmc)
#'
#' # Highlight clusters 2 and 3
#' getClusterColours(pbmc, clusters = c(2, 3))
#'
#' # Highlight clusters 2 and 3, set all other cluster colours to white
#' getClusterColours(pbmc, clusters = c(2, 3), default_colour = "white")
getClusterColours <- function(seurat, clusters = NULL,
                              original_colours = NULL, default_colour = "gray80") {

    # Handle clusters argument
    if (is.numeric(clusters)) clusters <- clusters + 1
    if (is.null(clusters)) clusters <- levels(seurat@ident)

    n_clust <- length(levels(seurat@ident))

    highlight_colours <- rep(default_colour, n_clust)
    names(highlight_colours) <- levels(seurat@ident)

    if (length(clusters) == 1) {
        if (clusters == "none") return(highlight_colours)
    }

    # Check if enough orig colours provided
    if (!is.null(original_colours) && length(original_colours) != length(levels(seurat@ident))) {

        if (length(original_colours) != length(clusters)) {

            stop("Not enough 'original_colours'! Please provide as many colours ",
                 "as clusters in the dataset, or one per cluster specified in the ",
                 "'clusters' argument.")

        } else highlight_colours[clusters] <- original_colours

    } else if (length(original_colours) == length(levels(seurat@ident))) {

        highlight_colours[clusters] <- original_colours[clusters]

    } else if (is.null(original_colours)) {

        original_colours <- ggColours(n_clust)
        names(original_colours) <- levels(seurat@ident)
        highlight_colours[clusters] <- original_colours[clusters]

    }

    return(highlight_colours)

}



#' Rotate the x axis labels in a ggplot
#'
#' @param angle Integer, value in degrees to rotate labels. Default: 90.
#'
#' @return A theme element to rotate labels
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + rotateX()
rotateX <- function(angle = 90) {

    theme(axis.text.x = element_text(angle = angle, hjust = 1))

}



#' Remove the legend in a ggplot
#'
#' @return A theme element to hide legend
#' @export
#'
#' @author Selin Jessa
#'
#' @examples
#' # gg <- mpg %>%
#' # filter(class %in% c("compact", "suv")) %>%
#' # ggplot(aes(x = displ, y = hwy)) +
#' # geom_point(aes(colour = factor(year))) +
#' # facet_wrap(~ class, ncol = 2)
#'
#' # gg
#' # gg + noLegend()
noLegend <- function() {

    theme(legend.position = "none")

}



#' Remove axis ticks and tick labels from a ggplot
#'
#' @return A theme element to remove ticks
#' @export
#'
#' @author Selin Jessa
noTicks <- function() {

    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

}


#' clusterCenters
#'
#' Get centers of clusters given a Seurat object, to use for labelling
#' in tSNE space. The cluster center is defined as the median X and Y coordinate
#' across cells in each cluster.
#'
#' @param seurat Seurat object, where dimensionality reduction has been applied,
#' i.e. (after applying Seurat::RunTSNE() to the object).
#'
#' @return Data frame with three columns: Cluster, mean_tSNE_1, and mean_tSNE_2
#' @export
#'
#' @author Selin Jessa
#' @examples
#'
#' clusterCenters(pbmc, reduction = "pca", dim1 = 1, dim2 = 3)
clusterCenters <- function(seurat, reduction, dim1, dim2) {

    n_clusters <- length(unique(seurat@ident))

    # Attempts at tidyeval...
    # vars <- colnames(seurat@dr[[reduction]]@cell.embeddings)[c(1, 2)]
    # col_names <- paste0("mean_", vars)

    # Get the embedding
    df <- as.data.frame(seurat@dr[[reduction]]@cell.embeddings[, c(dim1, dim2)]) %>%
        mutate(Cell = names(seurat@ident),
               Cluster = seurat@ident)

    # Generalize these
    colnames(df)[c(1, 2)] <- c("Dim_1", "Dim_2")

    # Compute cluster centers
    centers <- df %>%
        group_by(Cluster) %>%
        summarise(mean_x = median(Dim_1),
                  mean_y = median(Dim_2))

    return(centers)

}


#' Add cluster labels to a tSNE ggplot2 plot
#'
#' @param centers Data frame with at least three columns: "mean_x", "mean_y",
#' and "Cluster", as returned by \code{\link{clusterCenters}}
#' @param label_repel Logical, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#' @param label_short (Optional/Experimental!!) Logical, if TRUE, assumes clusters
#' (at \code{seurat@@ident}) consist of a prefix and a suffix separated by a non-alpha
#' numeric character (\code{"[^[:alnum:]]+"}), and tries to separate these names
#' and only plot the prefix, for shorter labels and a cleaner plot. Default: FALSE.
#' @param clusters (Optional) Clusters for which labels should be plot (if only
#' a subset of clusters should be labelled). Default: NULL (Label all clusters).
#'
#'
#' @author Selin Jessa
#' @export
addLabels <- function(centers, label_repel = FALSE, label_size = 4, label_short = FALSE, clusters = NULL) {

    if (!is.null(clusters)) centers <- filter(centers, Cluster %in% clusters)

    if (label_short) centers <- suppressWarnings(
        tidyr::separate(centers, Cluster, into = c("Cluster", "Cluster_long"), extra = "drop"))

    if (label_repel) {

        ggrepel::geom_label_repel(data = centers,
                                  aes(x = mean_x, y = mean_y),
                                  label = centers$Cluster,
                                  size = label_size,
                                  segment.color = 'grey50',
                                  fontface = 'bold',
                                  alpha = 0.8,
                                  segment.alpha = 0.8,
                                  label.size = NA,
                                  force = 2,
                                  # Leaving these unspecified for now, since it really depends on
                                  # the dimensionality reduction
                                  # nudge_x = 5, nudge_y = 5,
                                  segment.size = 0.5,
                                  arrow = arrow(length = unit(0.01, 'npc')))

    } else {

        geom_text(data = centers,
                  aes(x = mean_x, y = mean_y, label = Cluster),
                  size = label_size)

    }

}



#' Get the limits of a the first two dimensions in a dimensionality reduction
#'
#' When plotting an embedding, we may want to plot specific cells, but
#' constrain the scale to match plots of the whole dataset. Given a dim.
#' reduction, this function extracts the x and y limits to use for plotting.
#'
#' @param seurat Seurat object for which a dimensionality reduction has been
#' computed (e.g. PCA or tSNE)
#' @param reduction String, corresponding to the dimensionality reduction to use.
#' Default: "tsne".
#'
#' @return A list with two elements: "xlim", which is a character vector of
#' the limits for the x-axis, and "ylim", correspondingly for the y-axis
#' @export
#' @author Selin Jessa
#'
#' @examples
#' drLims(pbmc)
drLims <- function(seurat, reduction = "tsne") {

    dim1 <- seurat@dr[[reduction]]@cell.embeddings[,1]
    dim2 <- seurat@dr[[reduction]]@cell.embeddings[,2]

    return(list(xlim = c(min(dim1), max(dim1)),
                ylim = c(min(dim2), max(dim2))))

}




#' Constrain the scale of the plot to the dimensionality reduction limits
#'
#' @inheritParams drLims
#'
#' @export
#' @author Selin Jessa
constrainScale <- function(seurat, reduction = "tsne")  {

    limits <- drLims(seurat = seurat, reduction = reduction)
    lims(x = limits$xlim, y = limits$ylim)


}



#' Apply a clean theme to a ggplot2 object
#'
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
#' @author Selin Jessa
#' @export
theme_min <- function(base_size = 11, base_family = "",
                      border_colour = "grey90",
                      border_size = 1) {

    theme_light(base_size = 11, base_family = "") +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = border_colour, size = border_size),
            axis.ticks = element_line(colour = border_colour),
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



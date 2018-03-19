# Helper and utility functions for plotting



#' ggColours
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



#' rotateX
#'
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



#' noLegend
#'
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



#' noTicks
#'
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


#' addLabels
#'
#' Add cluster labels to a tSNE ggplot2 plot
#'
#' @param centers Data frame with at least three columns: "mean_tSNE_1", "mean_tSNE_2",
#' and "Cluster", as returned by \code{\link{clusterCenters}}
#' @param label_repel Logical, whether to plot cluster
#' labels repelled from the center, on a slightly transparent white background and
#' with an arrow pointing to the cluster center. If FALSE, simply plot the
#' cluster label at the cluster center. Default: TRUE.
#' @param label_size Numeric, controls the size of text labels. Default: 4.
#'
#' @author Selin Jessa
#' @export
addLabels <- function(centers, label_repel = FALSE, label_size = 4) {

    if (label_repel) {

        ggrepel::geom_label_repel(data = centers,
                                      aes(x = mean_tSNE_1, y = mean_tSNE_2),
                                      label = centers$Cluster,
                                      size = label_size,
                                      segment.color = 'grey50',
                                      fontface = 'bold',
                                      alpha = 0.8,
                                      segment.alpha = 0.8,
                                      label.size = NA,
                                      force = 2,
                                      nudge_x = 5, nudge_y = 5,
                                      segment.size = 0.5,
                                      arrow = arrow(length = unit(0.01, 'npc')))

    } else {

        geom_text(data = centers,
                  aes(x = mean_tSNE_1, y = mean_tSNE_2, label = Cluster),
                  size = label_size)

    }

}


#' theme_min
#'
#' A clean theme for ggplot2
#' @references https://github.com/sjessa/ggmin
#'
#' @importFrom ggplot2 theme_light theme
#' @author Selin Jessa
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

# Helper and utility functions

#' ggColours
#'
#' Get evenly spaced colours from around the colour wheel, which are the default
#' colours assigned to clusters by Seurat. The output of this function can be
#' passed to the \code{scale_colour_manual()} and \code{scale_fill_manual()} functions
#' from ggplot2, as the \code{values} argument.
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
ggColours <- function(n) {

    hues <- seq(15, 375, length = n + 1)
    colours <- hcl(h = hues, l = 65, c = 100)[1:n]
    names(colours) <- seq(0, n - 1) # Since the first cluster in Seurat is 0

    return(colours)

}

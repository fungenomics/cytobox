# Function to label tSNE clusters.

#' labelTSNEClusters
#'
#' @param seurat Seurat object.
#' @param labels A data frame with 3 columns: "cluster", "label" and "color". The "color" column is optional.
#' @param colors Cluster colors (Optional).
#'
#' @return A ggplot2 object. A tSNE plot with labelled clusters
#'
#' @export
#' @author Alexis Blanchet-Cohen
labelTSNEClusters <- function(seurat, labels) {

    old.cluster.ids <- levels(unique(GetClusters(seurat)$cluster))
    new.cluster.ids <- old.cluster.ids
    new.cluster.ids[labels$cluster+1] <- labels$label
    data@ident <- plyr::mapvalues(x = data@ident, from = old.cluster.ids, to = new.cluster.ids)

    if ("color" %in% labels$color) {
     # Selected colors.
        p <- TSNEPlot(object = seurat, do.label = FALSE, pt.size = 0.5, colors.use=labels$color, do.return=TRUE)
    } else {
        p <- TSNEPlot(object = seurat, do.label = FALSE, pt.size = 0.5, do.return=TRUE)
    }
    return(p)
}

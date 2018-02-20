# Function to label tSNE clusters.

#' labelTSNEClusters
#'
#' @param object Seurat object
#' @param labelsFile Path to labels file
#'
#' @return tSNE plot with labelled clusters
#'
#' @export
#' @author Alexis Blanchet-Cohen
labelTSNEClusters <- function(seurat, labelsFile) {

    labels <- data.table::fread(labelsFile)

    old.cluster.ids <- levels(unique(GetClusters(seurat)$cluster))
    new.cluster.ids <- old.cluster.ids
    new.cluster.ids[labels$cluster+1] <- labels$label

    if (colors %in% labels$color) {
     # Selected colors.
        p <- TSNEPlot(object = seurat, do.label = FALSE, pt.size = 0.5, colors.use=labels$color, do.return=TRUE)
    } else {
        p <- TSNEPlot(object = seurat, do.label = FALSE, pt.size = 0.5, do.return=TRUE)
    }
    return(p)
}

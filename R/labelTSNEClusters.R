# Plots a labelled tSNE plot from a Seurat object, where the clusters are labelled directly on the plot

#' labelTSNEClusters
#'
#' @param seurat Seurat object.
#' @param labels An data frame with 3 columns: "cluster", "label" and "colour". The "colour" column is optional. 
#' @param labelsOnPlot Boolean to indicate whether labels should be displayed on tSNE plot.
#' @param legend Boolean to indicate whether legend should be plotted.
#
#' @return A ggplot2 object. A tSNE plot with labelled clusters
#'
#' @export
#' @author Alexis Blanchet-Cohen
labelTSNEClusters <- function(seurat, labels = NULL, labelsOnPlot=TRUE, legend=TRUE) {
  # If labels are specified, replace old labels with new labels.
  if (!is.null(labels)) {
    old.cluster.ids <- levels(seurat@ident)
    new.cluster.ids <- old.cluster.ids
    new.cluster.ids[labels$cluster+1] <- labels$label
    seurat@ident <- plyr::mapvalues(x = seurat@ident, from = old.cluster.ids, to = new.cluster.ids)
  }

  # Get original Seurat colours
  n_clusters <- length(levels(seurat@ident))
  colours <- c(cytokit::ggColours(n_clusters))
  names(colours) <- levels(seurat@ident)

  # If colours are specified, replace Seurat colours with specified colours.
  if (!is.null(labels) & "colour" %in% colnames(labels)) {
    colours.dd <- data.frame(label=names(colours), colour=unname(colours))
    colours.dd <- left_join(labels, colours.dd, by="label")
    colours.dd <- mutate(colours.dd, colour.merged=ifelse(is.na(colour.x), colour.y, colour.x))
    colours <- colours.dd$colour.merged
    names(colours) <- colours.dd$label
  }

  # Get the embeddings
  df <- seurat@dr[["tsne"]]@cell.embeddings %>%
    as.data.frame %>%
    mutate(Cell = names(seurat@ident), Cluster = unname(seurat@ident)) %>%
    mutate(Cluster = factor(Cluster, levels = names(colours)))  

  # Compute cluster centers
  centers <- df %>%
    group_by(Cluster) %>%
    summarise(mean_tSNE_1 = mean(tSNE_1),
              mean_tSNE_2 = mean(tSNE_2))

  # Plot
  p1 <- df %>%
    ggplot(aes(x = tSNE_1, y = tSNE_2))

    p1 <- p1 +
      geom_point(aes(colour = Cluster), size = rel(0.8))

    if (label) {

      p1 <- p1 +
          ggrepel::geom_label_repel(data = centers,
                           aes(x = mean_tSNE_1, y = mean_tSNE_2, fill = Cluster),
                           label = centers$Cluster,
                           segment.colour = 'grey50',
                           force = 2,
                           nudge_x = 5, nudge_y = 5,
                           segment.size = 0.5,
                           arrow = arrow(length = unit(0.01, 'npc')))

    }

    p1 <- p1 +
      scale_colour_manual(values = colours) +
      scale_fill_manual(values = colours, guide = FALSE)

  p1 <- p1 + ggmin::theme_min()
  
  # Remove legend, if legend is set to FALSE.
  if(legend==FALSE) {
      p1 <- p1 + guides(colour=FALSE)
  }

  return(p1)

}

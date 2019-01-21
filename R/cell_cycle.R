#' computeCC
#'
#' Compute the cell cycle scores for a Seurat object, consisting of the mean expression,
#' in every cell, of a set of G1/S and G2/M phase markers.
#'
#' @param seurat Seurat object
#' @param species "mm" or "hg"
#'
#' @return A seurat object with scores in the \code{@meta.data} slot in two
#' new columns, "G1_S_score" and "G2_M_score"
#'
#' @export
#' @author Selin Jessa
computeCC <- function(seurat, species = "mm") {

    cc <- cytokit::cellCyclePlot(seurat, species = switch(species,
                                                          "mm" = "m_musculus",
                                                          "hg" = "h_sapiens"),
                                 return_scores = TRUE)

    seurat@meta.data$G1_S_score <- cc$g1.s.scores
    seurat@meta.data$G2_M_score <- cc$g2.m.scores

    return(seurat)


}

#' cellCyclePlot
#'
#' Plots a labelled tSNE plot from a Seurat object, where the clusters are labelled directly on the plot
#'
#' @param seurat Seurat object.
#' @param species m_musculus or h_sapiens.
#' @param facets Boolean to indicate whether a facet should be used for each cluster., or if all the clusters should be plotted together.
#' @param legend Boolean to indicate whether legend should be included.
#
#' @return A ggplot2 object. A tSNE plot with the cell cycle plots
#'
#' @export
#' @author Alexis Blanchet-Cohen
cellCyclePlot <- function(seurat, species="m_musculus", facets=TRUE, legend=FALSE,
                          return_scores = FALSE) {

  cell.cycle.genes <- cytobox::cell.cycle.genes.whitfield.2002
  if(species=="m_musculus") {
    cell.cycle.genes$gene.symbol <- cell.cycle.genes.whitfield.2002$mmusculus.gene.symbol
  } else {
    cell.cycle.genes$gene.symbol <- cell.cycle.genes.whitfield.2002$hsapiens.gene.symbol
  }
  g1.s.genes <- filter(cell.cycle.genes, phase=="G1/S") %>% .$gene.symbol
  g2.m.genes <- filter(cell.cycle.genes, phase=="G2/M") %>% .$gene.symbol

  expression.data <- as.data.frame(as.matrix(seurat@data))

  expression.data.g1.s.genes <- filter(expression.data, rownames(expression.data) %in% g1.s.genes)
  expression.data.g2.m.genes <- filter(expression.data, rownames(expression.data) %in% g2.m.genes)

  expression.data.g1.s.scores <- colMeans(expression.data.g1.s.genes)
  expression.data.g2.m.scores <- colMeans(expression.data.g2.m.genes)

  cell.cycle.scores <- as.data.frame(rbind(expression.data.g1.s.scores, expression.data.g2.m.scores))

  rownames(cell.cycle.scores) <- gsub("expression.data.", "", rownames(cell.cycle.scores))

  cell.cycle.scores.tidy <- as.data.frame(t(cell.cycle.scores))
  cell.cycle.scores.tidy <- tibble::rownames_to_column(cell.cycle.scores.tidy, "cell")
  cell.cycle.scores.tidy <- tibble::add_column(cell.cycle.scores.tidy, cluster=seurat@ident, .after="cell")

  if (return_scores) return(cell.cycle.scores.tidy)

  # Plots
  p <- ggplot(cell.cycle.scores.tidy, aes(x=g1.s.scores, y=g2.m.scores)) +
    geom_point(aes(color=cluster)) + xlab("G1/S score") + ylab("G2/M score")

  if(facets) {
  p <- p + facet_grid(~cluster) }

  if(legend==FALSE) {
    p <- p + theme(legend.position="none")
  }

  return(p)
}

# Plots a tSNE plot from a Seurat object, displaying the percentiles of the specified genes.

#' cellCyclePlot
#'
#' @param seurat Seurat object.
#' @param genes Genes to plot.
#' @param legend Boolean to indicate whether legend should be included.
#
#' @return A ggplot2 object. A tSNE plot with the cell cycle plots
#'
#' @export
#' @author Alexis Blanchet-Cohen
tsneByPercentileMarkerExpressionAlexis <- function(seurat, genes, legend=TRUE){
    # tSNE plots coordinates
    coordinates <- seurat@dr$tsne@cell.embeddings

    # Seurat normalized gene expression data
    expression.matrix <- as.data.frame(as.matrix(seurat@data))
    expression.matrix <- tibble::rownames_to_column(expression.matrix, "gene_name")
    # Keep only genes of interest.
    expression.matrix.filtered <- filter(expression.matrix, gene_name %in% genes)
    # Sum genes (No change if only one gene given).
    expression.matrix.filtered.summed <- colSums(select(expression.matrix.filtered, -gene_name))
    expression.matrix.filtered.summed <- data.frame(sum.genes=unname(expression.matrix.filtered.summed))

    coordinates.expression.data <- cbind(coordinates, expression.matrix.filtered.summed)

    # Add percentiles
    coordinates.expression.data.zero <- filter(coordinates.expression.data, sum.genes == 0) %>% mutate(percentile=NA)
    coordinates.expression.data.not.zero <- filter(coordinates.expression.data, sum.genes != 0)
    coordinates.expression.data.not.zero$percentile <- ecdf(coordinates.expression.data.not.zero$sum.genes)(coordinates.expression.data.not.zero$sum.genes) * 100
    coordinates.expression.data <- rbind(coordinates.expression.data.zero, coordinates.expression.data.not.zero)
    # Add color columns
    gradient.colors <- brewer.pal(n=9, name="Blues")
    coordinates.expression.data$gradient.colors.category <- 1
    coordinates.expression.data$gradient.colors.category <- ifelse(coordinates.expression.data$percentile > 0 & coordinates.expression.data$percentile <= 50,
						       2, coordinates.expression.data$gradient.colors.category)
    coordinates.expression.data$gradient.colors.category <- ifelse(coordinates.expression.data$percentile > 50 & coordinates.expression.data$percentile <= 70,
						       3, coordinates.expression.data$gradient.colors.category)
    coordinates.expression.data$gradient.colors.category <- ifelse(coordinates.expression.data$percentile > 70 & coordinates.expression.data$percentile <= 90,
						       4, coordinates.expression.data$gradient.colors.category)
    coordinates.expression.data$gradient.colors.category <- ifelse(coordinates.expression.data$percentile > 90 & coordinates.expression.data$percentile <= 92,
						       5, coordinates.expression.data$gradient.colors.category)
    coordinates.expression.data$gradient.colors.category <- ifelse(coordinates.expression.data$percentile > 92 & coordinates.expression.data$percentile <= 94,
						       6, coordinates.expression.data$gradient.colors.category)
    coordinates.expression.data$gradient.colors.category <- ifelse(coordinates.expression.data$percentile > 94 & coordinates.expression.data$percentile <= 96,
						       7, coordinates.expression.data$gradient.colors.category)
    coordinates.expression.data$gradient.colors.category <- ifelse(coordinates.expression.data$percentile > 96 & coordinates.expression.data$percentile <= 98,
						       8, coordinates.expression.data$gradient.colors.category)
    coordinates.expression.data$gradient.colors.category <- ifelse(coordinates.expression.data$percentile > 98 & coordinates.expression.data$percentile <= 100,
						       9, coordinates.expression.data$gradient.colors.category)
    coordinates.expression.data$gradient.colors.category <- ifelse(is.na(coordinates.expression.data$percentile),
							  1, coordinates.expression.data$gradient.colors.category)
      

    coordinates.expression.data$gradient.alpha.category <- 0.025
    coordinates.expression.data$gradient.alpha.category <- ifelse(coordinates.expression.data$percentile > 0 & coordinates.expression.data$percentile <= 50,
									  0.050, coordinates.expression.data$gradient.alpha.category)
    coordinates.expression.data$gradient.alpha.category <- ifelse(coordinates.expression.data$percentile > 50 & coordinates.expression.data$percentile <= 70,
									  0.10, coordinates.expression.data$gradient.alpha.category)
    coordinates.expression.data$gradient.alpha.category <- ifelse(coordinates.expression.data$percentile > 70 & coordinates.expression.data$percentile <= 90,
									  0.20, coordinates.expression.data$gradient.alpha.category)
    coordinates.expression.data$gradient.alpha.category <- ifelse(coordinates.expression.data$percentile > 90 & coordinates.expression.data$percentile <= 92,
									  0.30, coordinates.expression.data$gradient.alpha.category)
    coordinates.expression.data$gradient.alpha.category <- ifelse(coordinates.expression.data$percentile > 92 & coordinates.expression.data$percentile <= 94,
									  0.40, coordinates.expression.data$gradient.alpha.category)
    coordinates.expression.data$gradient.alpha.category <- ifelse(coordinates.expression.data$percentile > 94 & coordinates.expression.data$percentile <= 96,
									  0.50, coordinates.expression.data$gradient.alpha.category)
    coordinates.expression.data$gradient.alpha.category <- ifelse(coordinates.expression.data$percentile > 96 & coordinates.expression.data$percentile <= 98,
									  0.60, coordinates.expression.data$gradient.alpha.category)
    coordinates.expression.data$gradient.alpha.category <- ifelse(coordinates.expression.data$percentile > 98 & coordinates.expression.data$percentile <= 100,
									  1, coordinates.expression.data$gradient.alpha.category)
    coordinates.expression.data$gradient.alpha.category <- ifelse(is.na(coordinates.expression.data$percentile),
									 0.025, coordinates.expression.data$gradient.alpha.category)
     
    # Find values corresponding to percentiles.
    coordinates.expression.data <- dplyr::arrange(coordinates.expression.data, sum.genes)
    category.1.min <- category.1.max <- 0
    category.2.min <- filter(coordinates.expression.data, gradient.colors.category ==2) %>% slice(1) %>% .$sum.genes
    category.2.max <- filter(coordinates.expression.data, gradient.colors.category ==2) %>% slice(n()) %>% .$sum.genes
    category.3.min <- filter(coordinates.expression.data, gradient.colors.category ==3) %>% slice(1) %>% .$sum.genes
    category.3.max <- filter(coordinates.expression.data, gradient.colors.category ==3) %>% slice(n()) %>% .$sum.genes
    category.4.min <- filter(coordinates.expression.data, gradient.colors.category ==4) %>% slice(1) %>% .$sum.genes
    category.4.max <- filter(coordinates.expression.data, gradient.colors.category ==4) %>% slice(n()) %>% .$sum.genes
    category.5.min <- filter(coordinates.expression.data, gradient.colors.category ==5) %>% slice(1) %>% .$sum.genes
    category.5.max <- filter(coordinates.expression.data, gradient.colors.category ==5) %>% slice(n()) %>% .$sum.genes
    category.6.min <- filter(coordinates.expression.data, gradient.colors.category ==6) %>% slice(1) %>% .$sum.genes
    category.6.max <- filter(coordinates.expression.data, gradient.colors.category ==6) %>% slice(n()) %>% .$sum.genes
    category.7.min <- filter(coordinates.expression.data, gradient.colors.category ==7) %>% slice(1) %>% .$sum.genes
    category.7.max <- filter(coordinates.expression.data, gradient.colors.category ==7) %>% slice(n()) %>% .$sum.genes
    category.8.min <- filter(coordinates.expression.data, gradient.colors.category ==8) %>% slice(1) %>% .$sum.genes
    category.8.max <- filter(coordinates.expression.data, gradient.colors.category ==8) %>% slice(n()) %>% .$sum.genes
    category.9.min <- filter(coordinates.expression.data, gradient.colors.category ==9) %>% slice(1) %>% .$sum.genes
    category.9.max <- filter(coordinates.expression.data, gradient.colors.category ==9) %>% slice(n()) %>% .$sum.genes

    # Percentiles in legend
    cairo_pdf(file.path(outputDirectory, paste0(paste(genes, collapse="_"), "_distribution_umis_with_percentiles_in_legend.pdf")), width=10)
    p <- ggplot(coordinates.expression.data, aes(tSNE_1, tSNE_2)) +
        geom_point(aes(colour = factor(gradient.colors.category), alpha=gradient.alpha.category)) +
	scale_color_manual(values=gradient.colors, 
			   name="Expression level percentile\n(Color gradient)",
			   labels=c("Undetected", "> 0 & \u2264 50", "> 50 & \u2264 70", "> 70 & \u2264 90", "> 90 & \u2264 92", "> 92 & \u2264 94", "> 94 & \u2264 96", "> 96 & \u2264 98", "> 98 & \u2264 100")) +
    scale_alpha_continuous(name="Expression level percentile\n(Alpha gradient)", 
			       breaks=c(0.025, 0.050, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 1),
			       labels=c("Undetected", "> 0 & \u2264 50", "> 50 & \u2264 70", "> 70 & \u2264 90", "> 90 & \u2264 92", "> 92 & \u2264 94", "> 94 & \u2264 96", "> 96 & \u2264 98", "> 98 & \u2264 100")) +
   guides(colour = guide_legend(order = 1), 
        alpha = guide_legend(order = 2)) +
	xlab("tSNE_1") + ylab("tSNE_2") +
	ggtitle(paste(genes, collapse="_")) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5), legend.text.align = 0,
	      panel.grid.major = element_blank(), panel.grid.minor = element_blank())
      show(p)
      dev.off()

      # Expression values in legend
      if (legend == TRUE) {
      cairo_pdf(file.path(outputDirectory, paste0(paste(genes, collapse="_"), "_distribution_umis_with_expression_values_in_legend.pdf")), width=10)
      p <- ggplot(coordinates.expression.data, aes(tSNE_1, tSNE_2)) +
	geom_point(aes(colour = factor(gradient.colors.category), alpha=gradient.alpha.category)) +
	scale_color_manual(values=gradient.colors, 
			   name="Expression level\n(Color gradient)",
			   labels=c("Undetected", paste0("> 0 & \u2264 ", round(category.2.max, digits=2)), paste0("> ", round(category.2.max, digits=2), " & \u2264 ", round(category.3.max, digits=2)), paste0("> ", round(category.3.max, digits=2), " & \u2264 ", round(category.4.max, digits=2)), paste0("> ", round(category.4.max, digits=2), " & \u2264 ", round(category.5.max, digits=2)), paste0("> ", round(category.5.max, digits=2), " & \u2264 ", round(category.6.max, digits=2)), paste0("> ", round(category.6.max, digits=2), " & \u2264 ", round(category.7.max, digits=2)), paste0("> ", round(category.7.max, digits=2), " & \u2264 ", round(category.8.max, digits=2)), paste0("> ", round(category.8.max, digits=2), " & \u2264 ", round(category.9.max, digits=2)))) +
	scale_alpha_continuous(name="Expression level percentile\n(Alpha gradient)", 
			       breaks=c(0.025, 0.050, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 1),
			       labels=c("Undetected", paste0("> 0 & \u2264 ", round(category.2.max, digits=2)), paste0("> ", round(category.2.max, digits=2), " & \u2264 ", round(category.3.max, digits=2)), paste0("> ", round(category.3.max, digits=2), " & \u2264 ", round(category.4.max, digits=2)), paste0("> ", round(category.4.max, digits=2), " & \u2264 ", round(category.5.max, digits=2)), paste0("> ", round(category.5.max, digits=2), " & \u2264 ", round(category.6.max, digits=2)), paste0("> ", round(category.6.max, digits=2), " & \u2264 ", round(category.7.max, digits=2)), paste0("> ", round(category.7.max, digits=2), " & \u2264 ", round(category.8.max, digits=2)), paste0("> ", round(category.8.max, digits=2), " & \u2264 ", round(category.9.max, digits=2)))) +
	guides(colour = guide_legend(order = 1), 
	       alpha = guide_legend(order = 2)) +
	xlab("tSNE_1") + ylab("tSNE_2") +
	ggtitle(paste(genes,  collapse="_")) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5), legend.text.align = 0,
	      panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    } else {
      # No legend
      cairo_pdf(file.path(outputDirectory, paste0(paste(genes, collapse="_"), "_distribution_umis_with_expression_values_no_legend.pdf")), width=10)
      p <- ggplot(coordinates.expression.data, aes(tSNE_1, tSNE_2)) +
	geom_point(aes(colour = factor(gradient.colors.category), alpha=gradient.alpha.category)) +
	scale_color_manual(values=gradient.colors) +
	guides(color=FALSE, alpha=FALSE) +
	xlab("tSNE_1") + ylab("tSNE_2") +
	ggtitle(paste(genes, collapse="_")) +
	theme_bw() +
	theme(plot.title = element_text(hjust = 0.5), legend.text.align = 0,
	      panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    }
    return(p) 
}

# Monocle functions

#' Converts a Seurat object to a Monocle object.
#'
#' Full name of selected clusters must be given.
#' In your script, you must load the monocle library, before calling this function.
#'
#' @param seurat Seurat object to convert to a Monocle (CellDataSet) object.
#' @param selected_clusters Selected clusters in Seurat object. Full names of clusters must be given.
#'
#' @return A CellDataSet object
#' @export
#'
#' @author Alexis Blanchet-Cohen
#'
seurat_to_monocle <- function(seurat, selected_clusters) {
    # Parse seurat@meta.data to be in the correct format
    # Convert rownames with cell identity to column if it hasn't been done already.
    if (!"cell" %in% colnames(seurat@meta.data)) {
        seurat_meta_data <- rownames_to_column(seurat@meta.data, "cell")
        rownames(seurat_meta_data) <- seurat_meta_data$cell
        seurat@meta.data <- seurat_meta_data
    }
    # Copy res column to cluster column, if cluster column with labelling has not been created already.
    if (!"cluster" %in% colnames(seurat@meta.data)) {
        # Find first res column.
        index <- grep("res", colnames(seurat@meta.data))[1]
        seurat_meta_data$cluster <- seurat@meta.data[, index]
        rownames(seurat_meta_data) <- seurat_meta_data$cell
        seurat@meta.data <- seurat_meta_data
    }
    # Rename old column named percent.mito to percent_mito
    if ("percent.mito" %in% colnames(seurat@meta.data)) {
        seurat_meta_data <- rename(seurat@meta.data, percent_mito=percent.mito)
        rownames(seurat_meta_data) <- seurat_meta_data$cell
        seurat@meta.data <- seurat_meta_data
    }
    # Add default Seurat colours if the information is not in seurat@misc$colours
    if (is_empty(seurat@misc$colours)) {
        # Load the "scales" package
        suppressPackageStartupMessages(require(scales))
        # Create vector with levels of seurat@ident
        identities <- levels(seurat@ident)
        # Create vector of default ggplot2 colors
        my_color_palette <- hue_pal()(length(identities))
        seurat@misc$colours <- setNames(my_color_palette, identities)
    }
    # Keep only cells in selected clusters, and colours from selected clusters.
    if(!is.null(selected_clusters)){
        seurat_meta_data <- seurat@meta.data
        seurat_meta_data$cluster <- as.character(seurat_meta_data$cluster)
        if (!"cell" %in% colnames(seurat_meta_data)) { seurat_meta_data <- rownames_to_column(seurat_meta_data, "cell") }
        selected_cells <- filter(seurat_meta_data, cluster %in% selected_clusters) %>% .$cell
        seurat <- SubsetData(seurat, cells.use = selected_cells)
        seurat@meta.data$cluster <- as.factor(seurat@meta.data$cluster)
        seurat@meta.data$cluster <- droplevels(seurat@meta.data$cluster)
        seurat@misc$colours <- seurat@misc$colours[as.character(selected_clusters)]
    }
    cds <- monocle::importCDS(seurat)
    # Add the sample name.
    if (!is_empty(cds@experimentData@name)) { cds@experimentData@name <- seurat@project.name }
    # Add the misc data (colours and other information)
    cds@experimentData@other <- seurat@misc

    # Compute the size factors.
    cds <- BiocGenerics::estimateSizeFactors(cds)
    # Compute the dispersions. Slow.
    cds <- BiocGenerics::estimateDispersions(cds)
    # Detect the number of genes expressed.
    cds <- monocle::detectGenes(cds, min_expr = 0.1)

    # Computes a smooth function describing how variance in each gene's expression across cells varies according to the mean.
    disp_table <- monocle::dispersionTable(cds)

    # Ordering_genes
    ordering_genes <- as.character(base::subset(disp_table, mean_expression >= 0.1 &
                                 dispersion_empirical >= 1 * dispersion_fit)$gene_id)
    cds <- monocle::setOrderingFilter(cds, ordering_genes)

    # Reduce dimensions. Slow.
    cds <- monocle::reduceDimension(cds, residualModelFormulaStr="~num_genes_expressed + percent_mito")
    # Order cells
    cds <- monocle::orderCells(cds, reverse=FALSE)

    return(cds)
}


#' Draws a tSNE plot with a pseudotime gradient in the selected clusters.
#'
#'
#' @param seurat Seurat object to convert to a Monocle (CellDataSet) object.
#' @param monocle Monocle (CellDataSet) object with pseudotime data.
#' @param selected_clusters Selected clusters in Seurat object. Full names of clusters must be given.
#' @param legend. Boolean. Display legend or not.
#'
#' @return A ggplot2 object.
#' @export
#'
#' @author Alexis Blanchet-Cohen
#'
tsne_plot_clusters_highlighted_with_pseudotime <- function(seurat, monocle, selected_clusters, legend=TRUE){
    seurat_meta_data <- seurat@meta.data
    if(!"cell" %in% colnames(seurat_meta_data)) {
        seurat_meta_data <- tibble::rownames_to_column(seurat@meta.data, "cell")
        rownames(seurat_meta_data) <- seurat_meta_data$cell
        seurat@meta.data <- seurat_meta_data
    }
    # Copy first res column to cluster column, if cluster column with labelling has not been created already.
    if (!"cluster" %in% colnames(seurat@meta.data)) {
        # Find first res column.
        index <- grep("res", colnames(seurat@meta.data))[1]
        seurat_meta_data$cluster <- seurat@meta.data[, index]
        rownames(seurat_meta_data) <- seurat_meta_data$cell
        seurat@meta.data <- seurat_meta_data
    }
    seurat_meta_data$cluster <- as.character(seurat_meta_data$cluster)
    # Get the default colours if the colours column has not been set yet.
    if(is_empty(seurat@misc$colours)) {
        # Load the "scales" package
        suppressPackageStartupMessages(require(scales))
        # Create vector with levels of seurat@ident
        identities <- levels(seurat@ident)
        # Create vector of default ggplot2 colors
        my_color_palette <- hue_pal()(length(identities))
        seurat@misc$colours <- setNames(my_color_palette, identities)
    }
    seurat_meta_data <- seurat@meta.data

    monocle_PData <- Biobase::pData(monocle)
    if(!"cell" %in% colnames(monocle_PData)) { monocle_PData <- rownames_to_column(monocle_PData, "cell") }

    seurat_meta_data <- left_join(seurat_meta_data, dplyr::select(monocle_PData, cell, Pseudotime), by=c("cell"="cell"))

    rownames(seurat_meta_data) <- seurat_meta_data$cell
    seurat@meta.data <- seurat_meta_data

    # Personalized t-SNE plot.
    # tSNE plots coordinates
    coordinates <- seurat@dr$tsne@cell.embeddings
    coordinates <- cbind(coordinates, seurat@meta.data)
    coordinates$cluster <- as.character(coordinates$cluster)
    coordinates <- mutate(coordinates, cluster=ifelse(cluster %in% selected_clusters, cluster, "Not selected"))

    colours <- c(seurat@misc$colours[names(seurat@misc$colours) %in% selected_clusters], setNames("lightgrey", "Not selected"))

    alpha <- colours
    alpha[names(alpha) == "Not selected"] <- 0.1
    alpha[names(alpha) != "Not selected"] <- 1

    p <- ggplot(coordinates, aes(tSNE_1, tSNE_2, colour=Pseudotime, alpha=factor(cluster))) +
        geom_point() +
        scale_color_viridis(option="magma") +
        scale_alpha_manual(values=alpha) +
        guides(alpha=FALSE) +
        xlab("tSNE_1") + ylab("tSNE_2") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.text.align = 0,
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "right", legend.title=element_blank())

    if (legend==FALSE) {
        p <- p + guides(colour=FALSE, alpha=FALSE)
    }
    return(p)
}


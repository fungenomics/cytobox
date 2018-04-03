


# Random helpers ####

head2d <- function(mat, n = 10) {

    mat[1:n, 1:n]

}




# ---------------------------------------------------------------------------- #
# Re: clusters ####
# ---------------------------------------------------------------------------- #


whichCells <- function(seurat, clusters) {

    names(seurat@ident)[seurat@ident %in% clusters]

}

matchClusters <- function(seurat, pattern) {

    levels(seurat@ident)[which(str_detect(levels(seurat@ident), pattern))]

}

filterGenes <- function(seurat) {

    genes_hvg <- rownames(seurat@hvg.info)

    # Low expression: less than 20 detected molecules across cells
    n_umis_across_cells <- rowSums(as.matrix(seurat@raw.data))
    genes_low <- names(n_umis_across_cells[n_umis_across_cells < 20])

    genes_hvg_high <- setdiff(genes_hvg, genes_low)

    # Wide expression
    expr_bin <- as.matrix((as.matrix(seurat@data) > 0) + 0) # Converts numeric to binary matrix
    percent_detected <- rowSums(expr_bin)/ncol(expr_bin)
    genes_narrow <- names(percent_detected[percent_detected < 0.6])

    genes_hvg_high_narrow <- setdiff(genes_hvg_high, genes_narrow)

    return(genes_hvg_high_narrow)

}


clusterClusters <- function(seurat) {

    genes_filt <- filterGenes(seurat)
    exp_filt <- fetchData(seurat, genes_filt)
    exp_filt[exp_filt == 0] <- NA
    exp_means <- exp_filt %>%
        mutate(Cluster = seurat@ident) %>%
        group_by(Cluster) %>%
        summarise_all(funs(mean), na.rm = TRUE) %>%
        as.data.frame()

    # rownames(exp_means) <- as.character(rownames(exp_means))
    exp_means_noclust <- dplyr::select(exp_means, -Cluster)

    pheatmap::pheatmap(dist(exp_means_noclust),
                       # clustering_distance_rows = "correlation",
                       # clustering_distance_cols = "correlation",
                       clustering_method = "ward.D2",
                       labels_row = exp_means$Cluster,
                       labels_col = exp_means$Cluster,
                       color = viridis::magma(100),
                       dendrogram = "row")

}


relabelClusters <- function(seurat, labels) {

    old.cluster.ids <- levels(seurat@ident)
    seurat@ident <- plyr::mapvalues(x = seurat@ident, from = old.cluster.ids, to = labels)
    seurat@meta.data$cluster <- seurat@ident

    return(seurat)

}





subcluster <- function(seurat, clusters, genes, lab = TRUE) {

    # Get data
    df <- cytokit::fetchData(seurat, genes = genes, clusters = clusters,
                             return_cell = TRUE, return_cluster = TRUE,
                             scaled = FALSE)

    cluster_ids <- as.character(df$Cluster)
    rownames(df) <- df$Cell
    df <- dplyr::select(df, -Cell, -Cluster)
    mat <- t(as.matrix(df))

    if (!lab) {

        # Cluster
        gplots::heatmap.2(mat,
                          col = viridis::inferno(100),
                          labCol = FALSE,
                          trace = "none",
                          ylab = "Gene",
                          xlab = "Cell",
                          keysize = 1,
                          key.title = "Scale")


    } else {

        # Set colours
        colours <- ggColours(length(levels(seurat@ident)))
        cluster_colours <- cluster_ids

        for (i in seq_along(clusters)) {

            cluster_colours[cluster_colours == as.character(clusters[i])] <- colours[[as.character(clusters[i])]]

        }

        # Cluster
        gplots::heatmap.2(mat,
                          col = viridis::inferno(100),
                          labCol = FALSE,
                          trace = "none",
                          ColSideColors = cluster_colours,
                          ylab = "Gene",
                          xlab = "Cell",
                          keysize = 1,
                          key.title = "Scale")

    }



}





# ---------------------------------------------------------------------------- #
# Statistics and comparisons ####
# ---------------------------------------------------------------------------- #





# Following Hicks et al. (2017), CPM "normalizes each gene or transcript count by
# dividing the total number of UMIs per cell and multiplies by a scaling factor
# (10^6)
#
# @param obj: Seurat object
normalizeCPM <- function(obj) {

    mat <- as.matrix(obj@raw.data)
    (mat/colSums(mat)) * 10e6

}


# Calculate the cell-specific detection rate (adds as a column to the meta.data)
# @param obj: Seurat object
calcDetectionRate <- function(obj) {

    cell_x_gene <- t(normalizeCPM(obj))

    det_rate <- apply(cell_x_gene, 1, function(row) sum(row > 1)/ncol(cell_x_gene))

    keep_cells <- which((colnames(obj@raw.data) %in% rownames(obj@meta.data)))
    obj@meta.data$detection_rate <- det_rate[keep_cells]

    return(obj)

}



# !!! DEPRECATED !!! use cytokit::meanMarkerExprByCluster instead
#
# Adapted from Alexis' violin_plots.R script
# Given a set of clusters and their gene markers, compute the mean expression
# over all cells in a dataset for the markers for each cluster.
#
# @param clusters: Numeric vector of cluster numbers
# @param obj: Seurat object: does not have to be the same object from which the
# clusters/markers come from
# @param markers: Data frame of markers, generated by the pipeline, which should
# correspond to the clusters passed as the first argument
#
# @value: Data frame where the first column is the cell ID, the second is its
# cluster assignment in the obj dataset, and all others are
# the mean expression of marker genes for each cluster
meanMarkerExprByCluster <- function(clusters, obj, markers) {
    .Deprecated("cytokit::meanMarkerExprByCluster")

    perCluster <- function(cluster) {

        # Keep marker genes for cluster
        keep_markers <- markers[markers$cluster == cluster, ]$external_gene_name

        # Filtered expression
        filt_expr <- expr[rownames(expr) %in% keep_markers, ]

        # For each cell, compute its mean expression of gene markers for the cluster
        mean_expr <- colSums(filt_expr)/nrow(filt_expr)
        mean_expr <- data.frame(mean = as.numeric(mean_expr))

        names(mean_expr)[1] <- paste0(cluster)
        return(mean_expr)

    }

    expr <- obj@data %>% as.matrix %>% as.data.frame

    map_dfc(clusters, perCluster) %>%
        add_column(cluster = as.character(obj@ident), .before = 1) %>%
        add_column(cell = colnames(expr), .before = 1)

}


# !!! DEPRECATED !!! use cytokit::meanMarkerExpression instead
#
# Given a set of gene markers (say for a particular cell type, or cluster),
# compute each cell's mean expression over all the markers.
#
# @param obj: Seurat object
# @param markers: Character vector of gene names coresponding to the markers
#
# @value: Data frame with two columns: `cell`: the cell ID, and
# `mean_marker_expression` which contains the computed values
meanMarkerExp <- function(obj, markers) {
    .Deprecated("cytokit::meanMarkerExpression")

    # Expression data frome the Seurat object
    expr <- obj@data %>% as.matrix %>% as.data.frame

    # Filter to the markers
    filt_expr <- expr[rownames(expr) %in% markers, ]

    # Calculate mean expression
    mean_expr <- colSums(filt_expr)/nrow(filt_expr)
    mean_expr <- data.frame(cell = colnames(expr),
                            mean_marker_expression = as.numeric(mean_expr))

    return(mean_expr)

}

# !!! DEPRECATED !!! use cytokit::percentilesMarkerExpression instead
#
# Adapted from Alexis' scripts (`cell_type_markers_barreslab.R`)
calcExpressionPercentiles <- function(obj, markers) {
    .Deprecated("cytokit::percentilesMarkerExpression")

    expr.matrix <- obj@data %>%
        as.matrix() %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "gene_name")

    expr.matrix.filt <- filter(expr.matrix, gene_name %in% markers)
    expr.matrix.filt.summed <- colSums(select(expr.matrix.filt, -gene_name))
    expr.matrix.filt.summed <- data.frame(cell.type=unname(expr.matrix.filt.summed))

    # tSNE coordinates
    coord <- obj@dr$tsne@cell.embeddings
    coord.expr.data <- cbind(coord, expr.matrix.filt.summed)

    # Add percentiles
    coord.expr.data.zero <- filter(coord.expr.data, cell.type == 0) %>% mutate(percentile=NA)
    coord.expr.data.not.zero <- filter(coord.expr.data, cell.type != 0)
    coord.expr.data.not.zero$percentile <- ecdf(coord.expr.data.not.zero$cell.type)(coord.expr.data.not.zero$cell.type) * 100
    coord.expr.data <- rbind(coord.expr.data.zero, coord.expr.data.not.zero)

    # Add color columns
    # TODO try to replace this with a long case_when statement
    coord.expr.data$gradient.colors.group <- 1
    coord.expr.data$gradient.colors.group <- ifelse(coord.expr.data$percentile > 0 & coord.expr.data$percentile <= 50,
                                                    2, coord.expr.data$gradient.colors.group)
    coord.expr.data$gradient.colors.group <- ifelse(coord.expr.data$percentile > 50 & coord.expr.data$percentile <= 70,
                                                    3, coord.expr.data$gradient.colors.group)
    coord.expr.data$gradient.colors.group <- ifelse(coord.expr.data$percentile > 70 & coord.expr.data$percentile <= 90,
                                                    4, coord.expr.data$gradient.colors.group)
    coord.expr.data$gradient.colors.group <- ifelse(coord.expr.data$percentile > 90 & coord.expr.data$percentile <= 92,
                                                    5, coord.expr.data$gradient.colors.group)
    coord.expr.data$gradient.colors.group <- ifelse(coord.expr.data$percentile > 92 & coord.expr.data$percentile <= 94,
                                                    6, coord.expr.data$gradient.colors.group)
    coord.expr.data$gradient.colors.group <- ifelse(coord.expr.data$percentile > 94 & coord.expr.data$percentile <= 96,
                                                    7, coord.expr.data$gradient.colors.group)
    coord.expr.data$gradient.colors.group <- ifelse(coord.expr.data$percentile > 96 & coord.expr.data$percentile <= 98,
                                                    8, coord.expr.data$gradient.colors.group)
    coord.expr.data$gradient.colors.group <- ifelse(coord.expr.data$percentile > 98 & coord.expr.data$percentile <= 100,
                                                    9, coord.expr.data$gradient.colors.group)
    coord.expr.data$gradient.colors.group <- ifelse(is.na(coord.expr.data$percentile),
                                                    1, coord.expr.data$gradient.colors.group)

    coord.expr.data$gradient.alpha.group <- 0.025
    coord.expr.data$gradient.alpha.group <- ifelse(coord.expr.data$percentile > 0 & coord.expr.data$percentile <= 50,
                                                   0.050, coord.expr.data$gradient.alpha.group)
    coord.expr.data$gradient.alpha.group <- ifelse(coord.expr.data$percentile > 50 & coord.expr.data$percentile <= 70,
                                                   0.10, coord.expr.data$gradient.alpha.group)
    coord.expr.data$gradient.alpha.group <- ifelse(coord.expr.data$percentile > 70 & coord.expr.data$percentile <= 90,
                                                   0.20, coord.expr.data$gradient.alpha.group)
    coord.expr.data$gradient.alpha.group <- ifelse(coord.expr.data$percentile > 90 & coord.expr.data$percentile <= 92,
                                                   0.30, coord.expr.data$gradient.alpha.group)
    coord.expr.data$gradient.alpha.group <- ifelse(coord.expr.data$percentile > 92 & coord.expr.data$percentile <= 94,
                                                   0.40, coord.expr.data$gradient.alpha.group)
    coord.expr.data$gradient.alpha.group <- ifelse(coord.expr.data$percentile > 94 & coord.expr.data$percentile <= 96,
                                                   0.50, coord.expr.data$gradient.alpha.group)
    coord.expr.data$gradient.alpha.group <- ifelse(coord.expr.data$percentile > 96 & coord.expr.data$percentile <= 98,
                                                   0.60, coord.expr.data$gradient.alpha.group)
    coord.expr.data$gradient.alpha.group <- ifelse(coord.expr.data$percentile > 98 & coord.expr.data$percentile <= 100,
                                                   1, coord.expr.data$gradient.alpha.group)
    coord.expr.data$gradient.alpha.group <- ifelse(is.na(coord.expr.data$percentile),
                                                   0.025, coord.expr.data$gradient.alpha.group)

    # Find values corresponding to percentiles.
    coord.expr.data <- dplyr::arrange(coord.expr.data, cell.type)

    group.1.min <- group.1.max <- 0
    group.2.min <- filter(coord.expr.data, gradient.colors.group ==2) %>% dplyr::slice(1) %>% .$cell.type
    group.2.max <- filter(coord.expr.data, gradient.colors.group ==2) %>% dplyr::slice(n()) %>% .$cell.type
    group.3.min <- filter(coord.expr.data, gradient.colors.group ==3) %>% dplyr::slice(1) %>% .$cell.type
    group.3.max <- filter(coord.expr.data, gradient.colors.group ==3) %>% dplyr::slice(n()) %>% .$cell.type
    group.4.min <- filter(coord.expr.data, gradient.colors.group ==4) %>% dplyr::slice(1) %>% .$cell.type
    group.4.max <- filter(coord.expr.data, gradient.colors.group ==4) %>% dplyr::slice(n()) %>% .$cell.type
    group.5.min <- filter(coord.expr.data, gradient.colors.group ==5) %>% dplyr::slice(1) %>% .$cell.type
    group.5.max <- filter(coord.expr.data, gradient.colors.group ==5) %>% dplyr::slice(n()) %>% .$cell.type
    group.6.min <- filter(coord.expr.data, gradient.colors.group ==6) %>% dplyr::slice(1) %>% .$cell.type
    group.6.max <- filter(coord.expr.data, gradient.colors.group ==6) %>% dplyr::slice(n()) %>% .$cell.type
    group.7.min <- filter(coord.expr.data, gradient.colors.group ==7) %>% dplyr::slice(1) %>% .$cell.type
    group.7.max <- filter(coord.expr.data, gradient.colors.group ==7) %>% dplyr::slice(n()) %>% .$cell.type
    group.8.min <- filter(coord.expr.data, gradient.colors.group ==8) %>% dplyr::slice(1) %>% .$cell.type
    group.8.max <- filter(coord.expr.data, gradient.colors.group ==8) %>% dplyr::slice(n()) %>% .$cell.type
    group.9.min <- filter(coord.expr.data, gradient.colors.group ==9) %>% dplyr::slice(1) %>% .$cell.type
    group.9.max <- filter(coord.expr.data, gradient.colors.group ==9) %>% dplyr::slice(n()) %>% .$cell.type

    return(coord.expr.data)

}


# ---------------------------------------------------------------------------- #
#   Doing things with markers ####
# ---------------------------------------------------------------------------- #


intersectMarkers <- function(markers1, cluster1, markers2, cluster2) {

    m1 <- extractMarkers(markers1, cluster1, get = TRUE)
    m2 <- extractMarkers(markers2, cluster2, get = TRUE)
    base::intersect(m1, m2)

}


extractMarkers <- function(markers, cluster, get = FALSE, convert_fun = NULL) {

    i <- cluster
    mk <- markers %>% filter(cluster == i) %>% arrange(desc(avg_logFC))

    if (get) {

        if (is.null(convert_fun)) mk$external_gene_name
        else na.omit(convert_fun(mk$external_gene_name))
    }
    else mk

}



findMarkers <- function(pattern) {

    idx <- grep(pattern, names(markers_sj), ignore.case = TRUE)

    # Should return present and absent genes
    getGenes <- function(i) {

        present <- rownames(markers_sj)[which(markers_sj[[i]] == 1)]
        absent  <- rownames(markers_sj)[which(markers_sj[[i]] == -1)]

        list(present = present, absent = absent)

    }

    out <- lapply(idx, getGenes)
    names(out) <- names(markers_sj)[idx]

    return(out)

}





# ---------------------------------------------------------------------------- #
#   Manipulating Seurat objects ####
# ---------------------------------------------------------------------------- #

# !!! DEPRECATED !!! use cytokit::addEmbedding instead
#
# @param obj: Seurat object
# @param keep: Character vector of genes for which expression should be retrieved
#
# @value Data frame with tSNE coordinates (in columns tSNE_1 and tSNE_2),
# sample names, cell names, and expression for the given genes
getEmbedding <- function(obj, keep) {
    .Deprecated("cytokit::addEmbedding")


    df <- obj@dr$tsne@cell.embeddings %>%
        as.data.frame %>%
        mutate(sample = obj@meta.data$orig.ident,
               cell = rownames(obj@meta.data))

    expr <- obj@data[keep, ] %>%
        as.matrix %>%
        t() %>%
        as.data.frame

    df <- bind_cols(df, expr)

    return(df)

}

# !!! DEPRECATED !!! use cytokit::addEmbedding instead
#
# @param obj: Seurat object
#
# @value Data frame with tSNE coordinates (in columns tSNE_1 and tSNE_2),
# sample names, cell names, and expression for the given genes
getTSNEEmbedding <- function(obj) {
    .Deprecated("cytokit::clusterCenters")


    df <- obj@dr$tsne@cell.embeddings %>%
        as.data.frame %>%
        mutate(sample = obj@meta.data$orig.ident,
               cell = colnames(obj@data))

    return(df)

}


#' # !!! DEPRECATED !!! use cytokit::fetchData instead
#'
#' fetchData
#'
#' Similar to Seurat::FetchData except it doesn't thrown an error if a gene
#' is not found in the data. Also much more limited.
#'
#' @param seurat Seurat object
#' @param genes Genes to filter
#'
#' @return Expression matrix for genes specified
#'
#' @author Selin Jessa
fetchData <- function(seurat, genes) {
    .Deprecated("cytokit::fetchData")

    exp <- as.matrix(seurat@data)
    exp_filt <- as.data.frame(t(exp[which(rownames(exp) %in% genes),]))

    return(exp_filt)

}



# ---------------------------------------------------------------------------- #
#   Plotting ####
# ---------------------------------------------------------------------------- #

# violins ####

#' # !!! DEPRECATED !!! use cytokit::vln instead
#'
#' vln
#'
#' Similar to Seurat::VlnPlot() except it skips genes that are not found in the
#' data, and has an option to group plots by cluster instead of gene.
#'
#' @param seurat Seurat object
#' @param genes Genes to plot violins for
#' @param facet_by String, one of "gene" or "cluster". Default: "gene", clusters will
#' be on the x axis, and there will be one plot per gene (akin to Seurat::VlnPlot).
#' If "cluster", genes will be on the x axis, and there will be one plot per cluster.
#'
#' @return A ggplot2 object
#' @author Selin Jessa
vPlot <- function(seurat, genes, facet_by = "gene") {
    .Deprecated("cytokit::vln")

    exp <- as.matrix(seurat@data)

    # TODO This breaks if there is only one gene found from the list
    exp_filt <- as.data.frame(t(exp[which(rownames(exp) %in% genes),])) %>%
        tibble::add_column("Cluster" = seurat@ident, .before = 1) %>%
        gather(gene, expression, 2:length(.))

    if (facet_by == "gene") {

        gg <- ggplot(exp_filt, aes(x = Cluster, y = expression, fill = Cluster))

    } else if (facet_by == "cluster") {

        gg <- ggplot(exp_filt, aes(x = gene, y = expression, fill = Cluster))

    }

    gg <- gg +
        geom_violin(scale = "width", adjust = 1, trim = TRUE) +
        geom_jitter(size = 0.4, alpha = 0.5) +
        theme_min() +
        theme(legend.position = "none")

    if (facet_by == "gene") {

        gg <- gg + facet_wrap(~ gene, scales = "free_y", ncol = 2) +
            xlab("Gene") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1))

    } else if (facet_by == "cluster") {

        gg <- gg + facet_wrap(~ Cluster, scales = "free_y") +
            xlab("Identity")

    }

    return(gg)

}

VlnPlot <- function(seurat, genes,
                    from_sp = "mm",
                    to_sp = "mm",
                    ncol = ifelse(length(genes) %in% c(2, 4), 2, 3), rotate_x = FALSE) {

    Seurat::VlnPlot(seurat, cytokit::filterGenesSample(genes, seurat,
                                                       from_sp = from_sp, to_sp = to_sp),
                    point.size.use = 0.1, size.x.use = 0, nCol = ncol, x.lab.rot = rotate_x)

}



# heatmaps ####


# Plot a heatmap with clusters from two samples, where cells are coloured by each
# cluster's mean expression of markers defined by the other cluster
#
# @param s1, s2: Seurat objects, where data@ident stores the cluster assignments for cells
# @param markers1, markers2: Data frames as returned by Seurat::FindAllMarkers
# @param s1_name, s2_name: Strings giving sample names for s1 and s2, used in
# axes of the heatmap
#
# @value ggplot object, a heatmap using geom_raster
mutualMarkerExpHeatmap <- function(s1, s2, markers1, markers2, s1_name, s2_name) {

    s1_clusters <- unique(markers1$cluster)
    s2_clusters <- unique(markers2$cluster)

    message("Calculating median expression of markers in each cluster...")

    mdn_s1_exp_in_s2_clusters <- meanMarkerExprByCluster(s2_clusters, s1, markers2) %>%
        select(-cell) %>%
        group_by(cluster) %>%
        summarise_all(funs(median))

    mdn_s2_exp_in_s1_clusters <- meanMarkerExprByCluster(s1_clusters, s2, markers1) %>%
        select(-cell) %>%
        group_by(cluster) %>%
        summarise_all(funs(median))

    mdn_s1_exp_in_s2_clusters <- mdn_s1_exp_in_s2_clusters %>%
        mutate(cluster = as.numeric(cluster)) %>%
        arrange(cluster)

    message("Calculating minimum of the two matrices to find mutual high expression...")

    mdn_s2_exp_in_s1_clusters_t <- mdn_s2_exp_in_s1_clusters %>%
        as.data.frame %>%
        magrittr::set_rownames(.$cluster) %>%
        select(-cluster) %>%
        t() %>%
        as.data.frame %>%
        tibble::rownames_to_column(var = "cluster") %>%
        as.tibble()

    # Order columns properly
    tmp <- mdn_s2_exp_in_s1_clusters_t %>% select(-cluster)
    tmp <- tmp[order(as.numeric(names(tmp)))]
    mdn_s2_exp_in_s1_clusters_t <- tmp %>% tibble::add_column(mdn_s2_exp_in_s1_clusters_t$cluster, .before = 1)

    min_medians <- pmin(data.matrix(mdn_s1_exp_in_s2_clusters), data.matrix(mdn_s2_exp_in_s1_clusters_t)) %>%
        as.data.frame()

    message("Plotting heatmap...")

    gg <- min_medians %>% mutate(cluster = mdn_s1_exp_in_s2_clusters$cluster) %>%
        gather(s2_cluster, min_median_expr, 2:length(.)) %>%
        mutate(cluster = factor(cluster, levels = sort(s1_clusters, decreasing = TRUE)),
               s2_cluster = factor(s2_cluster, levels = sort(s2_clusters))) %>%
        ggplot(aes(x = s2_cluster, y = cluster)) +
        geom_raster(aes(fill = min_median_expr)) +
        geom_text(aes(label = round(min_median_expr, 2)), colour = "white", size = 3) +
        scale_x_discrete(position = "top") +
        scale_fill_gradientn(colors = viridis::viridis(100)) +
        theme_min() +
        xlab(s2_name) + ylab(s1_name)

    return(gg)

}


# tsne ####

# !!! DEPRECATED !!! use cytokit::tsneByPercentileMarkerExpression instead
#
# @param percentiles: Data frame as returned by calcExpressionPercentiles()
tsneByMarkerPercentiles <- function(percentiles, title, centers = NULL) {
    .Deprecated("cytokit::tsneByPercentileMarkerExpression")

    color_grad_labels <- c("Undetected", "> 0 & \u2264 50", "> 50 & \u2264 70", "> 70 & \u2264 90", "> 90 & \u2264 92", "> 92 & \u2264 94", "> 94 & \u2264 96", "> 96 & \u2264 98", "> 98 & \u2264 100")
    alpha_grad_labels <- color_grad_labels

    gg <- percentiles %>%
        ggplot(aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(colour = factor(gradient.colors.group)), alpha =  0.6) +
        #               alpha = gradient.alpha.group)
        scale_color_manual(#values = RColorBrewer::brewer.pal(n = 9, name = "Blues"),
            values = viridis(9),
            name = "Expression level percentile",
            labels = color_grad_labels) +
        # scale_alpha_continuous(name = "Expression level percentile\n(Alpha gradient)",
        #                        breaks = c(0.025, 0.050, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 1),
        #                        labels = alpha_grad_labels
        #  ) +
        guides(colour = guide_legend(order = 1),
               alpha = guide_legend(order = 2))

    if(!is.null(centers)) {

        gg <- gg + geom_label_repel(data = centers,
                                    aes(x = mean_tSNE_1, y = mean_tSNE_2),
                                    label = centers$Cluster,
                                    segment.color = 'grey50',
                                    force = 2,
                                    nudge_x = 5, nudge_y = 5,
                                    segment.size = 0.5,
                                    arrow = arrow(length = unit(0.01, 'npc')))
    }

    gg <- gg +
        xlab("tSNE_1") + ylab("tSNE_2") +
        ggtitle(title) +
        theme_min()

    return(gg)

}


# !!! DEPRECATED !!! use cytokit::tsneByMeanMarkerExpression instead
#
# Plot tSNE embedding coloured by mean expression of a group of marker genes
# (e.g. for a certain cell type or cluster)
#
# @param mean_exp: Data frame as returned by meanMarkerExp
# @param embedding: Data frame as returned by getEmbedding
tsneByMeanMarkerExp <- function(mean_exp, embedding, title, centers = NULL) {
    .Deprecated("cytokit::tsneByMeanMarkerExpression")

    df <- inner_join(mean_exp, embedding, by = "cell")

    gg <- df %>% ggplot(aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(colour = mean_marker_expression), size = rel(0.8)) +
        scale_color_viridis()

    if(!is.null(centers)) {

        gg <- gg + geom_label_repel(data = centers,
                                    aes(x = mean_tSNE_1, y = mean_tSNE_2),
                                    label = centers$Cluster,
                                    segment.color = 'grey50',
                                    force = 2,
                                    nudge_x = 5, nudge_y = 5,
                                    segment.size = 0.5,
                                    arrow = arrow(length = unit(0.01, 'npc')))
    }

    gg <- gg +
        theme_min() +
        ggtitle(title)

    return(gg)

}



# @param df1, df2: Embeddings as obtained from getEmbedding
# @param gene: String denoting the gene used to colour tSNE plots for each dataset
sideBySide <- function(df1, df2, gene) {

    p1 <- df1 %>%
        ggplot(aes(x = tSNE_1, tSNE_2)) +
        geom_point(aes_string(colour = gene), alpha = 0.5, size = rel(0.8)) +
        scale_color_viridis() +
        theme_min()

    p2 <- df2 %>%
        ggplot(aes(x = tSNE_1, tSNE_2)) +
        geom_point(aes_string(colour = gene), alpha = 0.5, size = rel(0.8)) +
        scale_color_viridis() +
        theme_min()

    p1 + p2

}




# !!! DEPRECATED !!! use cytokit::markerViolinPlot instead
#
# Make a grid of violin plots with the mean expression of each cluster, for markers
# of each cluster either from that same dataset, or from a different one
#
# @param expr: Data frame as returned by meanMarkerExprByCluster()
# @param sample1, sample2: Strings with sample IDs
markerViolinPlot <- function(expr, sample1, sample2, s1_clusters, s2_clusters) {
    .Deprecated("cytokit::pairwiseVln")

    clusters2 <- paste0(sample2, " cluster ", sort(s2_clusters))

    gg <- expr %>%
        gather(s2_cluster, mean_expression, 3:ncol(.)) %>%
        mutate(s1_cluster = factor(cluster, levels = sort(s1_clusters)),
               s2_cluster = paste0(sample2, " cluster ", s2_cluster)) %>%
        mutate(s2_cluster = factor(s2_cluster, levels = clusters2)) %>%
        ggplot(aes(x = s1_cluster, y = mean_expression)) +
        geom_violin(aes(fill = s2_cluster)) +
        facet_wrap(~ s2_cluster) +
        theme_min() +
        xlab(paste0(sample1, " cluster")) +
        scale_fill_discrete(name = paste0(sample2, " cluster")) +
        theme(panel.grid.major.x = element_line(colour = "grey90"))


    return(gg)

}


# !!! DEPRECATED !!! use cytokit::clusterCenters instead
#
# @value: Data frame with three columns: Cluster, mean_tSNE_1, mean_tSNE_2
# where the latter two can be used as x- and y-coordinates respectively
# to plot/label clusters at their centers
getTSNECenters <- function(obj, cluster_col) {
    .Deprecated("cytokit::clusterCenters")

    n_clusters <- length(unique(obj@meta.data[[cluster_col]]))

    # Get the embedding
    df <- obj@dr$tsne@cell.embeddings %>%
        as.data.frame %>%
        mutate(Cell = rownames(obj@meta.data),
               Cluster = obj@meta.data[[cluster_col]])

    if(stringr::str_detect(cluster_col, "res")) {

        df <- df %>%
            mutate(Cluster = factor(Cluster, levels = seq(0, n_clusters-1)))

    }


    # Compute cluster centers
    centers <- df %>%
        group_by(Cluster) %>%
        summarise(mean_tSNE_1 = mean(tSNE_1),
                  mean_tSNE_2 = mean(tSNE_2))

    return(centers)

}


getNewTSNECenters <- function(obj, cluster_col) {

    n_clusters <- length(unique(obj@meta.data[[cluster_col]]))

    # Get the embedding
    df <- obj@dr$tsne@cell.embeddings %>%
        as.data.frame %>%
        mutate(Cell = rownames(obj@meta.data),
               Cluster = obj@meta.data[[cluster_col]],
               Sample = obj@meta.data$orig.ident)

    if(stringr::str_detect(cluster_col, "res")) {

        df <- df %>%
            mutate(Cluster = factor(Cluster, levels = seq(0, n_clusters-1)))

    }


    # Compute cluster centers
    centers <- df %>%
        group_by(Sample, Cluster) %>%
        summarise(mean_tSNE_1 = mean(tSNE_1, trim = 0.15),
                  mean_tSNE_2 = mean(tSNE_2, trim = 0.15))

    return(centers)

}


# !!! DEPRECATED !!! use cytokit::tsne instead
#
# Plot a Labelled tSNE (lTSNE) from a Seurat object,
# where the clusters are labelled directly
# on the plot
#
# @param obj: Seurat object, with
# @param cluster_col: string giving the column in obj@meta.data containing the
# cluster assignments, e.g. "res.0.8"
# @param colour_by: string giving the column in obj@meta.data to colour points by.
# Default is to colour by cluster (the default value is "cluster")
# @param dr: Dimensionality reduction to use ("tsne" or "pca"), default: "tsne"
#
# @value ggplot object
lTSNE <- function(obj, cluster_col,
                  label = TRUE, label_repel = TRUE, colour_by = "cluster", dr = "tsne", group_colours = NULL) {
    .Deprecated("cytokit::tsne")

    # Get colours
    n_clusters <- length(unique(obj@meta.data[[cluster_col]]))
    colours <- c(gg_color_hue(n_clusters))
    names(colours) <- seq(0, n_clusters-1)

    # Get the embedding
    df <- obj@dr[[dr]]@cell.embeddings %>%
        as.data.frame %>%
        mutate(Cell = rownames(obj@meta.data),
               Cluster = obj@meta.data[[cluster_col]]) %>%
        mutate(Cluster = factor(Cluster, levels = names(colours)))

    if (colour_by != "cluster") {

        df[[colour_by]] <- obj@meta.data[[colour_by]]

    }

    # Compute cluster centers
    centers <- getTSNECenters(obj, cluster_col)

    # Plot
    p1 <- df %>%
        ggplot(aes(x = tSNE_1, y = tSNE_2))

    if (colour_by == "cluster") {

        p1 <- p1 +
            geom_point(aes(colour = Cluster), size = rel(0.8))

        if (label) {

            p1 <- p1 +
                geom_label_repel(data = centers,
                                 aes(x = mean_tSNE_1, y = mean_tSNE_2, fill = Cluster),
                                 label = centers$Cluster,
                                 segment.color = 'grey50',
                                 force = 2,
                                 nudge_x = 5, nudge_y = 5,
                                 segment.size = 0.5,
                                 arrow = arrow(length = unit(0.01, 'npc')))

        }

        p1 <- p1 +
            scale_colour_manual(values = colours) +
            scale_fill_manual(values = colours, guide = FALSE)

    } else {

        p1 <- p1 +
            geom_point(aes_string(colour = colour_by), size = rel(0.8))

        if (label) {

            p1 <- p1 +
                geom_label_repel(data = centers,
                                 aes(x = mean_tSNE_1, y = mean_tSNE_2),
                                 label = centers$Cluster,
                                 segment.color = 'grey50',
                                 force = 2,
                                 nudge_x = 5, nudge_y = 5,
                                 segment.size = 0.5,
                                 arrow = arrow(length = unit(0.01, 'npc')))

        }

        if(!is.null(group_colours)) {

            p1 <- p1 +
                scale_colour_manual(values = group_colours)

        }

    }

    p1 <- p1 + theme_min()
    return(p1)

}




# @param joint_df: Dataframe
# @param sample: string
# @param i: Cluster number or vector of cluster numbers to label
# @param n_clusters: integer, # of clusters in the sample
tsneHighlightCluster <- function(joint_df, sample, i, n_clusters) {
    .Deprecated("tsneHighlight")

    highlight_colours <- rep("gray80", n_clusters)
    names(highlight_colours) <- seq(0, n_clusters-1)

    orig_colours <- gg_color_hue(n_clusters)
    highlight_colours[i+1] <- orig_colours[i+1]

    p <- joint_df %>%
        filter(Sample == sample) %>%
        ggplot(aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes_string(fill = paste0(sample, "_cluster"),
                              colour = paste0(sample, "_cluster")), size = rel(0.8)) +
        theme_min() +
        scale_colour_manual(values = highlight_colours) +
        scale_fill_manual(values = highlight_colours) +
        ggtitle(paste0(sample, " cluster ", i))

    return(p)

}



# !!! DEPRECATED !!! use cytokit::tsne instead
#
tsne <- function(seurat, label = TRUE, point_size = 0.6, alpha = 0.8) {
    .Deprecated("cytokit::tsne")

    # Assuming that the order of the levels is correct in the seurat object,
    # this should find the colours of the original clusters, and whatever they've been renamed,
    # if and only if the number of new cluster IDs is equal to the number of old ones
    palette <- ggColors(length(levels(seurat@ident)))
    names(palette) <- levels(seurat@ident)
    centers <- clusterCenters(seurat)

    embedding <- data.frame(Cell = seurat@cell.names,
                            tSNE_1 = seurat@dr$tsne@cell.embeddings[, 1],
                            tSNE_2 = seurat@dr$tsne@cell.embeddings[, 2],
                            Cluster = seurat@ident,
                            stringsAsFactors = FALSE)

    gg <- ggplot(embedding, aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(colour = Cluster), size = point_size, alpha = alpha) +
        scale_color_manual(values = palette)

    if (label) {

        gg <- gg + ggrepel::geom_label_repel(data = centers,
                                             aes(x = mean_tSNE_1, y = mean_tSNE_2),
                                             label = centers$Cluster,
                                             segment.color = 'grey50',
                                             fontface = 'bold',
                                             alpha = 0.8,
                                             segment.alpha = 0.8,
                                             label.size = NA,
                                             force = 2,
                                             nudge_x = 5, nudge_y = 5,
                                             segment.size = 0.5,
                                             arrow = arrow(length = unit(0.01, 'npc')))

    }

    gg <- gg + theme_min()

    return(gg)


}




tsneHighlight <- function(seurat, clusters, original_colours = NULL, default_colour = "gray80", ...) {

    highlight_colours <- getClusterColours(seurat = seurat,
                                           clusters = clusters,
                                           original_colours = original_colours,
                                           default_colour = default_colour)

    cytokit::tsne(seurat, colours = highlight_colours, ...)

}

# tsneMerge <- function(seurat, genes, colours, point_size = 0.6) {
#
#     embedding <- cytokit::fetchData(seurat, genes) %>%
#         addEmbedding(seurat, exp_df, reduction = "tsne")
#
#     gg <- ggplot(embedding, aes(x = tSNE_1, y = tSNE_2))
#
#     for (i in seq_along(genes)) {
#
#         gg <- gg + geom_point(aes_string(colour = genes[i]), size = point_size, alpha = alpha) +
#
#
#     }
#
#
#
#
#
# }


addCellCycle <- function(seurat) {

    cc <- cytokit::cellCyclePlot(seurat, return_scores = TRUE)
    seurat@meta.data$G1_S_score <- cc$g1.s.scores
    seurat@meta.data$G2_M_score <- cc$g2.m.scores

    return(seurat)

}

#' Add the cell cycle scores to the seurat object
#'
#' Compute a score for G1-S and G2-M based on expression of genes associated
#' with these phases of the cell cycle, and add these as new columns of the
#' seurat object.
#'
#' @param seurat Seurat object
#' @param species Character, "m_musculus" or "h_sapiens"
#'
#' @return Seurat object, with two new columns in \code{seurat@meta.data}:
#' "G1_S_score" and "G2_M_score"
#' @export
#'
#' @author Selin Jessa
scoreCellCycle <- function(seurat, species = "m_musculus") {

    cc <- cellCyclePlot(seurat, species = species, return_scores = TRUE)
    seurat@meta.data$G1_S_score <- cc$g1.s.scores
    seurat@meta.data$G2_M_score <- cc$g2.m.scores

    return(seurat)

}




# helpers ####

ggRotateX <- function(angle = 90) {
    .Deprecated("cytokit::rotateX")

    theme(axis.text.x = element_text(angle = angle, hjust = 1))

}

ggNoLegend <- function() {
    .Deprecated("cytokit::noLegend")

    theme(legend.position = "none")

}


cowplotAddTitle <- function(plots, title, ncol = 3) {

    plot_list <- cowplot::plot_grid(plotlist = plots, ncol = ncol)

    plot_title <- cowplot::ggdraw() + cowplot::draw_label(title, hjust = 1, size = 12)
    cowplot::plot_grid(plot_title, plot_list, ncol = 1, rel_heights = c(0.05, 1))

}

getClusterColours <- function(seurat, clusters, original_colours = NULL, default_colour = "gray80") {

    n_clust <- length(levels(seurat@ident))

    highlight_colours <- rep(default_colour, n_clust)
    names(highlight_colours) <- levels(seurat@ident)

    if (is.null(original_colours)) original_colours <- cytokit::ggColours(n_clust)
    names(original_colours) <- levels(seurat@ident)
    highlight_colours[clusters] <- original_colours[clusters]

    return(highlight_colours)

}



# !!! DEPRECATED !!! use cytokit::ggColours instead
#
# From Alexis
gg_color_hue <- function(n) {
    .Deprecated("cytokit::ggColours")

    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}





# Multi-sample plotting ####


#' @param joint_seurat Seurat object
#' @param indiv_seurat Named list of seurat objects
jointDf <- function(joint_seurat, indiv_seurat) {

    getCluster <- function(cell, sample_i) {

        raw_cell <- str_sub(cell, -16, -1)
        varhandle::unfactor(indiv_seurat[[sample_i]]@ident[[raw_cell]])

    }

    # Make a new column which contains the original clusters
    df <- data.frame(Cell = joint_seurat@cell.names,
                     Sample = joint_seurat@meta.data$orig.ident,
                     stringsAsFactors = FALSE)

    samples <- unique(joint_seurat@meta.data$orig.ident)
    df$Original_cluster <- NA

    # For each sample, get original cluster for all cells, and save it
    # to a new column
    for (i in 1:length(samples)) {

        idx_in_sample <- which(df$Sample == samples[i])

        df[[samples[i]]] <- NA

        df[idx_in_sample, ][[samples[i]]] <- lapply(df$Cell[idx_in_sample],
                                                    getCluster,
                                                    samples[i])

        df[[samples[i]]] <- factor(df[[samples[i]]],
                                   levels = levels(indiv_seurat[[samples[i]]]@ident))

        # Also create one column which contains each cell's sample of origin
        # and original cluster assignment
        df[idx_in_sample, ][["Original_cluster"]] <- paste0(df[idx_in_sample, ]$Sample,
                                                            "_",
                                                            df[idx_in_sample, ][[samples[i]]])

    }

    df$Sample <- factor(df$Sample, levels = samples)
    return(df)

}


# A general function to take in a joint Seurat object
# and generate some useful plots to visualize the clusters
tsneJoint <- function(joint_seurat, indiv_seurat,
                      point_size = 0.5, ncol = 2, trim = 0.15, alpha = 0.8,
                      label = TRUE, label_short = FALSE, return_list = FALSE, legend = TRUE) {

    # Get the dataframe for plotting
    df <- jointDf(joint_seurat, indiv_seurat)
    samples <- unique(joint_seurat@meta.data$orig.ident)

    # Get the joint embedding
    df <- addEmbedding(joint_seurat, df, reduction = "tsne")

    # Get new centers
    centers <- df %>% group_by(Sample, Original_cluster) %>%
        summarise(mean_tSNE_1 = mean(tSNE_1, trim = trim),
                  mean_tSNE_2 = mean(tSNE_2, trim = trim))
    centers$Cluster <- str_split(centers$Original_cluster, "_") %>%
        purrr::map(tail, 1) %>%
        unlist

    tsneSample <- function(sample_i) {

        p1 <- filter(df, !is.na(!!rlang::sym(sample_i))) %>%
            ggplot(aes(x = tSNE_1, y = tSNE_2)) +
            geom_point(aes_string(colour = sample_i, fill = sample_i), size = point_size, alpha = alpha) +
            theme_min() +
            ggtitle(sample_i)

        centers_sample <- filter(centers, Sample == sample_i)
        p1 <- p1 + cytokit::addLabels(centers_sample, label_repel = TRUE, label_short = label_short)

        if (!legend) p1 <- p1 + cytokit::noLegend()
        return(p1)

    }

    plot_list <- lapply(samples, tsneSample)

    if (return_list) return(plot_list)
    else cowplot::plot_grid(plotlist = plot_list, ncol = ncol)

}




tsneHighlightJoint <- function(joint_seurat, indiv_seurat, sample,
                               clusters, default_colour = "lightsteelblue1", label_short = FALSE) {

    seurat <- indiv_seurat[[sample]]
    n_clust <- length(levels(seurat@ident))

    highlight_colours <- rep(default_colour, n_clust)
    names(highlight_colours) <- levels(seurat@ident)

    orig_colours <- cytokit::ggColours(n_clust)
    names(orig_colours) <- levels(seurat@ident)
    highlight_colours[clusters] <- orig_colours[clusters]

    # Get the dataframe for plotting
    df <- jointDf(joint_seurat, indiv_seurat)
    samples <- unique(joint_seurat@meta.data$orig.ident)

    # Get the joint embedding
    df <- addEmbedding(joint_seurat, df, reduction = "tsne")

    # Get new centers
    centers <- df %>% group_by(Sample, Original_cluster) %>%
        summarise(mean_tSNE_1 = mean(tSNE_1, trim = 0.15),
                  mean_tSNE_2 = mean(tSNE_2, trim = 0.15))

    centers$Cluster <- str_split(centers$Original_cluster, "_") %>%
        purrr::map(tail, 1) %>%
        unlist

    p1 <- filter(df, !is.na(!!rlang::sym(sample))) %>%
        ggplot(aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes_string(colour = sample, fill = sample), size = 0.5) +
        scale_color_manual(values = highlight_colours) +
        theme_min() +
        ggtitle(sample) +
        cytokit::noLegend()

    centers_sample <- filter(centers, Sample == sample)
    p1 <- p1 + cytokit::addLabels(centers_sample, label_repel = TRUE, label_short = label_short)

    return(p1)

}


scale <- function(x) {
    (x - min(x)) / (max(x) - min(x))
}


#' usage: rgbGradient(ct_e15_proj, c("Pax6", "Eomes", "Rbfox3"))
rgbGradient <- function(seurat, genes, clusters = NULL, point_size = 2, alpha = 0.6) {

    genes <- head(genes, 3)

    # Get data
    df <- cytokit::fetchData(seurat, genes, clusters, return_cell = TRUE) %>%
        arrange_(genes[1]) %>%
        # Scale and convert to RGB
        mutate(R = scale(.[[2]]),
               G = scale(.[[3]]),
               B = scale(.[[4]])) %>%
        # Get hex code from RGB
        mutate(rgb = rgb(R, G, B))

    colours <- df$rgb
    names(colours) <- df$Cell

    cytokit::addEmbedding(seurat, df, reduction = "tsne") %>%
        ggplot(aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(aes(colour = Cell), stroke = 0, size = point_size, alpha = alpha) +
        scale_colour_manual(values = colours) +
        theme_min() + xlab("tSNE 1") + ylab("tSNE 2") +
        cytokit::noLegend()

}




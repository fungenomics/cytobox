# Functions for pairwise comparisons of two SC datasets



#' meanMarkerExprByCluster
#'
#' Compute the mean expression of markers for each cluster in one dataset,
#' in each cluster of another dataset. Used for \code{\link{markerViolinPlot}}
#'
#' @return Data frame where the first two columns are "Cell" and "Cluster",
#' and all remaining columns give the mean expression of markers identified
#' for each cluster
#'
#' @export
#' @author Selin Jessa
meanMarkerExprByCluster <- function(seurat, markers, marker_col = "gene") {

    perCluster <- function(i) {

        # Keep marker genes for cluster
        keep_markers <- filter(markers, cluster == i) %>% `[[`(marker_col)

        # Filtered expression
        filt_exp <- exp[rownames(exp) %in% keep_markers, ]

        # For each cell, compute its mean expression of gene markers for the cluster
        mean_exp <- colSums(filt_exp)/nrow(filt_exp)
        mean_exp <- data.frame(mean = as.numeric(mean_exp))

        names(mean_exp)[1] <- i
        return(mean_exp)

    }

    exp <- as.data.frame(as.matrix(seurat@data))
    clusters <- sort(unique(markers$cluster))

    # Bind columns, each of which corresponds to a cluster
    purrr::map_dfc(clusters, perCluster) %>%
        tibble::add_column(Cluster = as.character(seurat@ident), .before = 1) %>%
        tibble::add_column(Cell = colnames(exp), .before = 1)

}


#' markerViolinPlot
#'
#' Make a grid of violin plots with the mean expression of each cluster, for markers
# of each cluster either from that same dataset, or from a different one. Here,
#' dataset 1 is the sample whose expression values will be plot, and dataset 2
#' is the sample whose cluster markers are used.
#'
#' @param seurat1 Seurat object for which to plot expression
#' @param markers Data frame as returned by Seurat::FindAllMarkers()
#' @param seurat2 Seurat object from which cluster markers (\code{markers}) were defined
#' @param s1_name String, sample name for \code{seurat1}
#' @param s2_name String, sample name for \code{seurat2}
#' @param marker_col String specifying the column in the \code{markers} data frames which
#' specifies the cluster. By default, Seurat calls this "gene" (the default here);
#' in objects produced by the lab's pipeline, it may be called "external_gene_name".
#'
#' @return A ggplot object
#' @author adapted from Alexis Blanchet-Cohen
#'
#' @export
#' @examples
#' markerViolinPlot(pbmc, markers_pbmc, pbmc, "Test1", "Test2")
markerViolinPlot <- function(seurat1, markers, seurat2, s1_name, s2_name, marker_col = "gene") {

    s1_clusters <- sort(unique(seurat1@ident))
    s2_clusters <- sort(unique(seurat2@ident))
    clusters2 <- paste0(s2_name, " cluster ", s2_clusters)

    exp <- meanMarkerExprByCluster(seurat1, markers, marker_col)

    gg <- exp %>%
        tidyr::gather(s2_cluster, mean_expression, 3:ncol(.)) %>%
        dplyr::mutate(s1_cluster = factor(Cluster, levels = s1_clusters),
               s2_cluster = glue("{s2_name} cluster {s2_cluster}")) %>%
        mutate(s2_cluster = factor(s2_cluster, levels = clusters2)) %>%
        ggplot(aes(x = s1_cluster, y = mean_expression)) +
        geom_violin(aes(fill = s2_cluster)) +
        facet_wrap(~ s2_cluster) +
        theme_min() +
        xlab(glue("{s1_name} cluster")) +
        scale_fill_discrete(name = glue("{s2_name} cluster")) +
        theme(panel.grid.major.x = element_line(colour = "grey90"))

    return(gg)

}



#' percentMarkerOverlap
#'
#' @param markers1 Data frame of markers for one dataset as returned by Seurat::FindAllMarkers()
#' @param markers2 Data frame of markers for second dataset as returned by Seurat::FindAllMarkers()
#' @param marker_col String specifying the column in the markers data frames which
#' specifies the cluster. By default, Seurat calls this "gene" (the default here);
#' in objects produced by the lab's pipeline, it may be called "external_gene_name".
#'
#' @return Data frame with three columns: "s1_cluster" giving the cluster for
#' dataset 1, "s2_cluster" giving the cluster for dataset 2, and "marker_overlap"
#' which gives the computed min percent overlap of markers between the two clusters
#' @export
#' @author Selin Jessa
#' @examples
#' percentMarkerOverlap(markers_pbmc, markers_pbmc)
percentMarkerOverlap <- function(markers1, markers2, marker_col = "gene") {

    n_clust1 <- length(unique(markers1$cluster))
    n_clust2 <- length(unique(markers2$cluster))

    olaps = vector("list", n_clust1 * n_clust2)
    k <- 1

    for (i in 0:(n_clust1 - 1)) {

        for (j in 0:(n_clust2 - 1)) {

            n_olap <- length(intersect(markers1[markers1$cluster == i, marker_col],
                                markers2[markers2$cluster == j, marker_col]))

            olaps[[k]] <- min(n_olap / length(markers1[markers1$cluster == i, marker_col]),
                              n_olap / length(markers2[markers2$cluster == j, marker_col]))

            k <- k + 1

        }
    }

    min_olap <- data.frame(s1_cluster = rep(0:(n_clust1-1), each = n_clust2),
                           s2_cluster = rep(0:(n_clust2-1), times = n_clust1),
                           marker_overlap = unlist(olaps))

    return(min_olap)

}



#' heatmapPercentMarkerOverlap
#'
#' Generate a heatmap where rows correspond to clusters from one dataset and
#' columns correspond to the clusters from the second dataset. The value at [i, j]
#' is the minimum of the percent overlap between markers for Cluster i in dataset 1
#' and Cluster j in dataset 2.
#'
#' @param markers1 Data frame of markers for one dataset as returned by Seurat::FindAllMarkers()
#' @param markers2 Data frame of markers for second dataset as returned by Seurat::FindAllMarkers()
#' @param s1_name (Optional) string giving the name or ID of sample 1 (used for axis labels)
#' @param s2_name (Optional) string giving the name or ID of sample 2 (used for axis labels)
#' @param marker_col String specifying the column in the markers data frames which
#' specifies the cluster. By default, Seurat calls this "gene" (the default here); in the pipeline,
#' it may be called "external_gene_name".
#'
#' @export
#' @author Selin Jessa
#'
#' @examples
#' # Here, we compute the heatmap for one dataset (pbmc) with itself
#' heatmapPercentMarkerOverlap(markers_pbmc, markers_pbmc)
heatmapPercentMarkerOverlap <- function(markers1, markers2,
                                        s1_name = NULL, s2_name = NULL,
                                        marker_col = "gene") {

    min_olap <- percentMarkerOverlap(markers1, markers2)

    # Get nice ordering
    olap_srt <- min_olap %>% arrange(desc(marker_overlap))

    gg <- min_olap %>%
        dplyr::mutate(cluster = factor(s1_cluster, levels = rev(unique(olap_srt$s1_cluster))),
               s2_cluster = factor(s2_cluster, levels = unique(olap_srt$s2_cluster))) %>%
        ggplot(aes(x = s2_cluster, y = cluster)) +
        geom_raster(aes(fill = marker_overlap)) +
        scale_fill_gradientn(colors = viridis::viridis(100)) +
        geom_text(aes(label = round(marker_overlap, 2)), colour = "white", size = 3) +
        scale_x_discrete(position = "top") +
        theme_min() +
        theme(panel.border = element_blank()) +
        ggtitle("Min % overlap in cluster gene markers, sorted by % overlap")

    if ((!is.null(s1_name)) & (!is.null(s2_name))) {

        gg <- gg + xlab(s2_name) + ylab(s1_name)

    }

    return(gg)

}



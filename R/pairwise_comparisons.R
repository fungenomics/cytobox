# Functions for pairwise comparisons of two SC datasets


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
        geom_text(aes(label = round(marker_overlap, 2)), colour = "white", size = 3) +
        scale_x_discrete(position = "top") +
        scale_fill_viridis_c(limits = c(0, 1)) +
        theme_min() +
        theme(panel.border = element_blank()) +
        ggtitle("Min % overlap in cluster gene markers, sorted by % overlap")

    if ((!is.null(s1_name)) & (!is.null(s2_name))) {

        gg <- gg + xlab(s2_name) + ylab(s1_name)

    }

    return(gg)

}



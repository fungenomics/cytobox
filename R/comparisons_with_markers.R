# Functions for pairwise comparisons of two SC datasets



#' meanClusterMarkerExpr
#'
#' For each cluster, compute each cell's mean expression of that cluster's
#' markers. Used for \code{\link{markerViolinPlot}}.
#'
#' @param seurat Seurat object, whose expression values will be used
#' @param markers Markers data frame, for the same or another Seurat object. The
#' expression of markers in each cluster in this data frame will be calculated
#' for each cluster in \code{seurat}.
#'
#' @return Data frame where the first two columns are "Cell" and "Cluster" (from
#' \code{seurat}),
#' and all remaining columns give the mean expression of markers identified
#' for each cluster in \code{markers}
#'
#' @export
#' @author Selin Jessa
#' @examples
#'
#' library(dplyr)
#'
#' # Using the same sample's marker:
#' meanClusterMarkerExpr(pbmc, markers_pbmc, "gene")
#'
#' # Change the name of the clusters in the markers df, as if it were
#' # from a different sample where the clusters are A, B, C, D:
#' markers2 <- mutate(markers_pbmc, cluster = recode(
#'     cluster, `0` = "A", `1` = "B", `2` = "C", `3` = "D"))
#'
#' meanClusterMarkerExpr(pbmc, markers2, "gene")
meanClusterMarkerExpr <- function(seurat, markers, marker_col = "gene") {

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
    # Iterate over the clusters for which we have markers
    clusters <- as.character(sort(unique(markers$cluster)))

    # Bind columns, each of which corresponds to a cluster
    purrr::map_dfc(clusters, perCluster) %>%
        # Create a column with the clusters from the Seurat object
        tibble::add_column(Cluster = as.character(seurat@ident), .before = 1) %>%
        tibble::add_column(Cell = colnames(exp), .before = 1)

}


#' pairwiseVln
#'
#' Make a grid of violin plots with the mean expression of each cluster, for markers
# of each cluster either from that same dataset, or from a different one. Here,
#' dataset 1 is the sample whose expression values will be plot, and dataset 2
#' is the sample whose cluster markers are used.
#'
#' @param seurat1 Seurat object for which to plot expression
#' @param markers Data frame as returned by Seurat::FindAllMarkers()
#' @param seurat2 Seurat object from which cluster markers (\code{markers}) were defined
#' @param sample_names (Optional) Character vector giving the names or ID of
#' the two samples, used for axis labels, etc.
#' Default: c(seurat1@@project.name, seurat2@@project.name)
#' @param marker_col String specifying the column in the \code{markers} data frames which
#' specifies the cluster. By default, Seurat calls this "gene" (the default here);
#' in objects produced by the lab's pipeline, it may be called "external_gene_name".
#'
#' @return A ggplot object
#' @author adapted from Alexis Blanchet-Cohen
#'
#' @aliases markerViolinPlot
#' @export
#' @examples
#' pairwiseVln(pbmc, markers_pbmc, pbmc)
pairwiseVln <- function(seurat1, markers, seurat2,
                        sample_names = c(seurat1@project.name, seurat2@project.name),
                        marker_col = "gene") {

    clusters2 <- glue("{sample_names[2]} cluster {levels(seurat2@ident)}")

    exp <- meanClusterMarkerExpr(seurat1, markers, marker_col)

    gg <- exp %>%
        tidyr::gather(s2_cluster, mean_expression, 3:ncol(.)) %>%
        dplyr::mutate(s1_cluster = factor(Cluster, levels = levels(seurat1@ident)),
               s2_cluster = glue("{sample_names[2]} cluster {s2_cluster}")) %>%
        mutate(s2_cluster = factor(s2_cluster, levels = clusters2)) %>%
        ggplot(aes(x = s1_cluster, y = mean_expression)) +
        geom_violin(aes(fill = s2_cluster)) +
        facet_wrap(~ s2_cluster) +
        theme_min() +
        xlab(glue("{sample_names[1]} cluster")) +
        scale_fill_discrete(name = glue("{sample_names[2]} cluster")) +
        theme(panel.grid.major.x = element_line(colour = "grey90"))

    return(gg)

}


#' @export
markerViolinPlot <- pairwiseVln



#' percentMarkerOverlap
#'
#' @param markers1 Data frame of markers for one dataset as returned by Seurat::FindAllMarkers()
#' @param markers2 Data frame of markers for second dataset as returned by Seurat::FindAllMarkers()
#' @param mode "s1", "s2", "min", "max"
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
percentMarkerOverlap <- function(markers1, markers2, mode = "min", marker_col = "gene") {

    n_clust1 <- length(unique(markers1$cluster))
    n_clust2 <- length(unique(markers2$cluster))

    olaps_s1 <- vector("list", n_clust1 * n_clust2)
    olaps_s2 <- vector("list", n_clust1 * n_clust2)
    k <- 1

    for (i in 0:(n_clust1 - 1)) {

        for (j in 0:(n_clust2 - 1)) {

            mk_1 <- markers1[markers1$cluster == i,][[marker_col]]
            mk_2 <- markers2[markers2$cluster == j,][[marker_col]]
            n_olap <- length(base::intersect(mk_1, mk_2))
            olaps_s1[[k]] <- n_olap / length(mk_1)
            olaps_s2[[k]] <- n_olap / length(mk_2)
            k <- k + 1

        }
    }

    df_olap <- data.frame(s1_cluster = rep(0:(n_clust1-1), each = n_clust2),
                           s2_cluster = rep(0:(n_clust2-1), times = n_clust1))

    if (mode == "s1") df_olap$marker_overlap <- unlist(olaps_s1)
    else if (mode == "s2") df_olap$marker_overlap <- unlist(olaps_s2)
    else if (mode == "min") df_olap$marker_overlap <- pmin(unlist(olaps_s1), unlist(olaps_s2))
    else if (mode == "max") df_olap$marker_overlap <- pmax(unlist(olaps_s1), unlist(olaps_s2))

    return(df_olap)

}



#' Compute a heatmap of marker overlap between samples
#'
#' Generate a heatmap where rows correspond to clusters from one dataset and
#' columns correspond to the clusters from the second dataset. The value at [\emph{i, j}]
#' is overlap of markers between Cluster \emph{i} in dataset 1 and Cluster \emph{j} in dataset 2.
#' See Details for how to control how the "overlap" is computed.
#'
#' @details
#' When calculating the shared markers between clusters in two samples, we typically
#' count the number of markers in both lists, but have two options for what to use
#' as the denominator (either the number of markers in sample 1, or the number of
#' markers in sample 2). This function provides a few options for what to plot in the heatmap,
#' controlled by the `mode` argument:
#'
#' \itemize{
#'   \item "s1": Plot the percentage of sample 1 markers which are also markers in sample 2
#'   \item "s2": Plot the percentage of sample 2 markers which are also markers in sample 1
#'   \item "min": Plot the minimum of the above percentages
#'   \item "max": Plot the maximum of the above percentages
#'   \item "both": Plot two heatmaps, one with the percentage of sample 1 markers which are also markers in sample 2,
#' and one the percentage of sample 2 markers which are also markers in sample 1
#' }
#'
#' @param markers1 Data frame of markers for one dataset as returned by Seurat::FindAllMarkers()
#' The clusters for this sample will be plot on the x axis.
#' @param markers2 Data frame of markers for second dataset as returned by Seurat::FindAllMarkers()
#' The clusters for this sample will be plot on the y axis.
#' @param mode String, one of "s1", "s2", "min", "max", or "both". Controls how to
#' compute the marker overlaps. Default: "min". See Details for more info.
#' @param sample_names (Optional) Character vector giving the names or ID of
#' the two samples, used for axis labels, etc. Default: c("Sample 1", "Sample 2")
#' @param marker_col String specifying the column in the markers data frames which
#' specifies the cluster. By default, Seurat calls this "gene"; in the pipeline,
#' it may be called "external_gene_name" (the default here).
#' @param palette Character vector containing a gradient palette to use.
#' Default: \code{\link{viridis::magma}}.
#' @param label_colour String specifying the colour of value of each cell
#' printed in the heatmap. Default: "white"
#'
#' @export
#' @author Selin Jessa
#'
#' @examples
#' # Compute the heatmap for one dataset (pbmc) with itself
#' heatmapPercentMarkerOverlap(markers_pbmc, markers_pbmc, marker_col = "gene")
heatmapPercentMarkerOverlap <- function(markers1, markers2, mode = "min",
                                        sample_names = c("Sample 1", "Sample 2"),
                                        palette = NULL,
                                        marker_col = "external_gene_name",
                                        label_colour = "white") {

    if (mode == "both") {

        return(plot_grid(heatmapPercentMarkerOverlap(markers1, markers2, "s1",
                                                     sample_names = sample_names,
                                              palette = palette, marker_col = marker_col),
                  heatmapPercentMarkerOverlap(markers1, markers2, "s2",
                                              sample_names = sample_names,
                                              palette = palette, marker_col = marker_col),
                  ncol = 2))

    }

    min_olap <- percentMarkerOverlap(markers1, markers2, mode = mode, marker_col = marker_col)

    # Get nice ordering
    olap_srt <- min_olap %>% arrange(desc(marker_overlap))

    if (is.null(palette)) pal <- viridis::magma(100)
    else pal <- palette

    gg <- min_olap %>%
        dplyr::mutate(cluster = factor(s1_cluster, levels = rev(unique(olap_srt$s1_cluster))),
               s2_cluster = factor(s2_cluster, levels = unique(olap_srt$s2_cluster))) %>%
        ggplot(aes(x = s2_cluster, y = cluster)) +
        geom_raster(aes(fill = marker_overlap)) +
        scale_fill_gradientn(colors = pal, limits = c(0, 1)) +
        geom_text(aes(label = round(marker_overlap, 2)), colour = label_colour, size = 3) +
        scale_x_discrete(position = "top") +
        theme_min() +
        theme(panel.border = element_blank()) +
        xlab(sample_names[1]) + ylab(sample_names[2]) +
        guides(fill = guide_legend(title = "Marker overlap")) +
        ggtitle(case_when(
            mode == "min" ~ "Minimum % marker overlap",
            mode == "max" ~ "Maximum % marker overlap",
            mode == "s1" ~ as.character(glue("% {sample_names[1]} markers in {sample_names[2]}")),
            mode == "s2" ~ as.character(glue("% {sample_names[2]} markers in {sample_names[1]}"))))

    return(gg)

}



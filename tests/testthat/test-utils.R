context("test-utils.R")

test_that("getting ggplot2 default colours works", {
    cols <- c("0" = "#F8766D", "1" = "#00BFC4")
    expect_equal(ggColours(2), cols)
})


test_that("we can flexibly add embedding to a dataframe", {

    df <- data.frame(Cell = rownames(pbmc@meta.data),
                                     Cluster = pbmc@meta.data$res.0.8,
                     stringsAsFactors = FALSE)

    tsne_emb <- addEmbedding(pbmc, df)
    expect_equal(names(tsne_emb), c("Cell", "Cluster", "tSNE_1", "tSNE_2"))

    pca_emb <- addEmbedding(pbmc, df, reduction = "pca")
    expect_equal(names(pca_emb), c("Cell", "Cluster", "PC1", "PC2"))

})

test_that("addEmbedding matches data per cell correctly", {

    df <- data.frame(Cell = rownames(pbmc@meta.data),
                     Cluster = pbmc@meta.data$res.0.8,
                     stringsAsFactors = FALSE)

    tsne_emb <- addEmbedding(pbmc, df)

    expect_equal(tsne_emb[1, c("tSNE_1", "tSNE_2")],
                 data.frame(tSNE_1 = 14.42191, tSNE_2 = 8.336022),
                 tolerance = 0.1)

})


test_that("subsetting and fetching expression data works", {

    expect_error(fetchData(pbmc, "foo"), "No genes specified")
    expect_equal(nrow(fetchData(pbmc, c("IL32"), c(1, 2))), 41)
    expect_equal(nrow(fetchData(pbmc, c("IL32", "MS4A1"))), 80)
    expect_equal(ncol(fetchData(pbmc, c("IL32", "MS4A1"))), 2)
    expect_equal(names(fetchData(pbmc, c("IL32"), return_cluster = TRUE, return_cell = TRUE)),
                 c("Cell", "Cluster", "IL32"))

})


test_that("findGenes works", {

    find_out <- findGenes(pbmc, c("IL32", "CD79B", "foo"))
    expect_equal(find_out$detected, c("IL32", "CD79B"))
    expect_equal(find_out$undetected, "foo")

})



# test_that("we get cluster centers in tSNE space correctly", {
#
#     df <- tibble(Cluster = factor(c("0", "1", "2", "3")),
#                  mean_tSNE_1 = c(-1.14, 6.51, -17.1, 5.03),
#                  mean_tSNE_2 = c(19.6, 0.512, -5.65, -51.6))
#
#     expect_equal(clusterCenters(pbmc), df, tolerance = 1e-2)
#
# })

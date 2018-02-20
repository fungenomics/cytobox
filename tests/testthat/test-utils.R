context("test-utils.R")

test_that("getting ggplot2 default colours works", {
    cols <- c("0" = "#F8766D", "1" = "#00BFC4")
    expect_equal(ggColours(2), cols)
})


test_that("we can flexibly add embedding to a dataframe", {

    df <- data.frame(Cell = rownames(pbmc@meta.data),
                                     Cluster = pbmc@meta.data$res.0.8)

    tsne_emb <- addEmbedding(pbmc, df)
    expect_equal(names(tsne_emb), c("Cell", "Cluster", "tSNE_1", "tSNE_2"))

    pca_emb <- addEmbedding(pbmc, df, reduction = "pca")
    expect_equal(names(pca_emb), c("Cell", "Cluster", "PC1", "PC2"))

})

test_that("addEmbedding matches data per cell correctly", {

    df <- data.frame(Cell = rownames(pbmc@meta.data),
                     Cluster = pbmc@meta.data$res.0.8)

    tsne_emb <- addEmbedding(pbmc, df)

    expect_equal(tsne_emb[1, c("tSNE_1", "tSNE_2")],
                 data.frame(tSNE_1 = 14.42191, tSNE_2 = 8.336022),
                 tolerance = 1e-6)

})

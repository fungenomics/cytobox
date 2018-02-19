context("test-summarise.R")

test_that("meanMarkerExpression works", {

    one_gene <- data.frame(
        Cell = factor(c("ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC"),
                      levels = sort(rownames(pbmc@meta.data))),
        Mean_marker_expression = c(4.968821, 0, 6.942623)
    )

    expect_equal(head(meanMarkerExpression(pbmc, "IL32"), 3), one_gene,
                 tolerance = 1e-6)

    two_genes <- data.frame(
        Cell = factor(c("ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC"),
                      levels = sort(rownames(pbmc@meta.data))),
        Mean_marker_expression = c(4.968821, 0, 6.192272)
    )

    expect_equal(head(meanMarkerExpression(pbmc, c("IL32", "CD2")), 3), two_genes,
                 tolerance = 1e-6)

})

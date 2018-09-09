context("test-summarise.R")

test_that("meanGeneExpression works", {

    one_gene <- data.frame(
        Cell = c("ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC"),
        Mean_marker_expression = c(4.968821, 0.000000, 6.942623), stringsAsFactors = FALSE)

    expect_equal(head(meanGeneExpression(pbmc, "IL32"), 3), one_gene,
                 tolerance = 1e-6)

    two_genes <- data.frame(
        Cell = c("ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC"),
        Mean_marker_expression = c(4.968821, 0.000000, 6.192272), stringsAsFactors = FALSE)

    expect_equal(head(meanGeneExpression(pbmc, c("IL32", "CD2")), 3), two_genes,
                 tolerance = 1e-6)

})

context("test-correlations.R")


test_that("meanClusterExpression works", {

    mat <- matrix(c(2.563624, 1.239862, 0.5996828, 0.3811947,
                    2.6203414, 0.7697423, 0.3782335, 0.0000000), nrow = 2)
    rownames(mat) <- c("IL32", "CD2")
    colnames(mat) <- levels(pbmc@ident)

    expect_equal(meanClusterExpression(pbmc, c("IL32", "CD2")), mat,
                 tolerance = 1e-6)

})



test_that("correlations are computed correctly", {

    out <- matrix(
        c(1.0000000, 0.5192243, 0.6113777, 0.3836883,
          0.5192243, 1.0000000, 0.4611809, 0.4211603,
          0.6113777, 0.4611809, 1.0000000, 0.2909938,
          0.3836883, 0.4211603, 0.2909938, 1.0000000),
        nrow = 4)

    rownames(out) <- colnames(out) <- levels(pbmc@ident)

    expect_equal(correlateExpression(pbmc, pbmc, head(rownames(pbmc@raw.data), 100), "hg", "hg"),
                 out, tolerance = 1e-6)

})

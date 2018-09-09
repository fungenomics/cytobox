context("test-comparisons.R")

test_that("meanMarkerExpressionPerCell works", {

    markers2 <- dplyr::mutate(markers_pbmc,
                              cluster = recode(cluster,
                                               `0` = "A", `1` = "B", `2` = "C", `3` = "D"))

    df_out <- data.frame(Cell = c("ATGCCAGAACGACT", "CATGGCCTGTGCAT", "GAACCTGATGAACC"),
                         Cluster = c("0", "0", "0"),
                         A = as.numeric(c(1.962009, 2.362410, 2.482422)),
                         B = as.numeric(c(1.015172, 1.256989, 1.270847)),
                         C = as.numeric(c(0.9229822, 1.3528365, 0.5587361)),
                         D = as.numeric(c(1.488204, 1.086767, 1.013582)),
                         stringsAsFactors = FALSE)

    expect_equal(head(meanMarkerExpressionPerCell(pbmc, markers2), 3),
                 df_out, tolerance = 1e-3)

})

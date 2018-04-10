context("test-plotting_utils.R")

test_that("getting cluster colours works", {
   expect_equal(getClusterColours(pbmc, clusters = c(2, 3)),
                c("0" = "gray80", "1" = "gray80", "2" = "#00BFC4", "3" = "#C77CFF"))

    expect_error(getClusterColours(pbmc, clusters = c(2, 3), original_colours = "black"),
                 "Not enough")


})


test_that("retrieval of correct dimensionality reduction limits", {

    lims <- list(xlim = c(-21.29, 19.82),
                 ylim = c(-55.73, 33.26))

    expect_equal(drLims(pbmc), lims, tolerance = 1e-2)

})

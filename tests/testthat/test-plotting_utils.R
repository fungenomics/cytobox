context("test-plotting_utils.R")

test_that("getting cluster colours works", {
   expect_equal(getClusterColours(pbmc, clusters = c(2, 3)),
                c("0" = "gray80", "1" = "gray80", "2" = "#00BFC4", "3" = "#C77CFF"))

    expect_error(getClusterColours(pbmc, clusters = c(2, 3), original_colours = "black"),
                 "Not enough")


})

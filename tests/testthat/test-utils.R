context("test-utils.R")

test_that("getting ggplot2 default colours works", {
    cols <- c("0" = "#F8766D", "1" = "#00BFC4")
    expect_equal(ggColours(2), cols)
})

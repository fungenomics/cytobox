context("test-genes.R")

test_that("gene conversion works", {

    expect_equal(hg2mm("TUBB", "NES"), c("Nes", "Tubb5"))
    expect_equal(hg2mm(c("TUBB", "NES")), c("Nes", "Tubb5"))

})

test_that("gene filtering works", {
    expect_equal(filterGenesSample(c("IL32", "foo"), pbmc, "hg", "hg"), "IL32")
})

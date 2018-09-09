context("test-genes.R")

test_that("gene conversion works", {

    expect_equal(hg2mm("TUBB", "NES"), c("Nes", "Tubb5"))
    expect_equal(hg2mm(c("TUBB", "NES")), c("Nes", "Tubb5"))

})

test_that("gene filtering works", {
    expect_equal(filterGenesSample(c("IL32", "foo"), pbmc, "hg", "hg"), "IL32")
})


test_that("conversion between gene symbols and ENSEMBL IDs", {

    # 1A. symbols2ensembl on a matrix
    mat1 <- matrix(seq(1, 15), nrow = 3)
    rownames(mat1) <- c("Gfap", "Fabp7", "Foo")

    mat2 <- mat1[1:2, ]
    rownames(mat2) <- c("ENSMUSG00000020932", "ENSMUSG00000019874")

    expect_equal(symbols2ensembl(mat1, sp = "mm"), mat2)

    # 1B. on a vector
    expect_equal(symbols2ensembl(c("Gfap", "Fabp7", "SDKFJDF"), "mm"),
                 c("ENSMUSG00000020932", "ENSMUSG00000019874"))

    # 2. ensembl2symbols
    expect_equal(ensembl2symbols(c("ENSMUSG00000020932", "ENSMUSG00000019874"), "mm"),
                 c("Gfap", "Fabp7"))

})

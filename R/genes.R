# Functions for dealing with gene names
# Written Claudia Kleinman, adapted by Selin Jessa


#' geneName
#'
#' Find the proper capitalization (which is species dependent) for genes
#' TODO: return object does not preserve order of original list
#'
#' @param symbols String, or character vector of gene names
#' @param species Species for which \code{symbols} are given
#'
#' @return Rows of annotation corresponding to the symbols
#' @export
#'
#' @examples
#' geneName("TUBB", species = "hg")
geneName <- function(symbols, species = c("hg", "mm")) {
    switch(species, "hg" = hg.genes, "mm" = mm.genes) %>%
        filter(toupper(gene_symbol) %in% toupper(symbols))
}



#' convertName
#'
#' @param symbol String or character vector of gene names
#' @param from One of "hg" or "mm"
#' @param to  One of "hg" or "mm"
#'
#' @return Character vector of converted gene names
#' @export
#'
#' @examples
#' convertName(c("TUBB", "NES"), "hg", "mm")
convertName <- function(symbol, from, to = NULL) {
    switch(from, "hg" = hg2mm.genes, "mm" = mm2hg.genes) %>%
        filter(toupper(gene_symbol) %in% toupper(symbol)) %>%
        .$homologous_gene_symbol
}


#' mm2hg
#'
#' Wrapper around \code{\link{convertName}}
#'
#' @param ... Mouse genes to convert, as individual strings or one character vector
#'
#' @return Character vector of converted names
#' @export
#'
#' @examples
#' mm2hg("Tubb5", "Nes")
mm2hg <- function(...) {
    convertName(c(...), "mm")
}


#' hg2mm
#'
#' Wrapper around \code{\link{convertName}}
#'
#' @param ... Human genes to convert, as individual strings or one character vector
#'
#' @return Character vector of converted names
#' @export
#'
#' @examples
#' hg2mm("TUBB", "NES")
#' hg2mm(c("TUBB", "NES"))
hg2mm <- function(...) {
    convertName(c(...), "hg")
}


#' filterGenesSample
#'
#' @param genes Genes to search for
#' @param seurat Seurat object
#' @param from_sp One of "hg" or "mm". Default: "hg"
#' @param to_sp One of "hg" or "mm". Default: "mm"
#'
#' @return Character vector of converted gene names which are detected in the sample
#' @export
#'
#' @examples
#' filterGenesSample(c("IL32", "foo"), pbmc, "hg", "hg")
filterGenesSample <- function(genes, seurat, from_sp = "hg", to_sp = "mm") {

    sample_genes <- rownames(seurat@hvg.info)
    if (from_sp == to_sp) {
        filter(geneName(genes, from_sp), gene_symbol %in% sample_genes) %>% .$gene_symbol
    }
    else {
        conv_genes <- convertName(genes, from_sp)
        conv_genes[which(conv_genes %in% sample_genes)]
    }

}

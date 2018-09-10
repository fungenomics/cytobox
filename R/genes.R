# Functions for dealing with gene names
# Claudia Kleinman, Selin Jessa


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
#' @author adapted from Claudia Kleinman
#'
#' @examples
#' geneName("TUBB", species = "hg")
geneName <- function(symbols, species = c("hg", "mm")) {
    switch(species, "hg" = hg.genes, "mm" = mm.genes) %>%
        filter(toupper(gene_symbol) %in% toupper(symbols))
}

#' convertGenesBetweenSpeciesHomologs
#'
#' @param x Character vector or matrix with mouse gene symbols to convert.
#' If a matrix, rownames should be gene symbols.
#' @param from One of "hg" or "mm"
#' @param to  One of "hg" or "mm"
#'
#' @author adapted from Claudia Kleinman
#' @export
#'
#' @examples
#' convertGenesBetweenSpeciesHomologs(c("TUBB", "NES"), "hg", "mm")
#'
#' mat1 <- matrix(seq(1, 15), nrow = 3)
#' rownames(mat1) <- c("Gfap", "Fabp7", "Foo")
#' convertGenesBetweenSpeciesHomologs(mat1, from = "mm")
convertGenesBetweenSpeciesHomologs <- function (x, ...) {
    UseMethod("convertGenesBetweenSpeciesHomologs", x)
}

#' mm2hg
#'
#' Convert genes from mouse to human, wrapper around \code{\link{convertGenesBetweenSpeciesHomologs}}.
#' For vectors, the output is \strong{not order preserving}.
#' For matrices, order is preserved, but \strong{genes where the homolog is not
#' found are silently dropped}.
#'
#' @param x Character vector or matrix with mouse gene symbols to convert.
#' If character, provide either individual strings or one character vector.
#' If a matrix, rownames should be gene symbols.
#'
#' @export
#' @author Selin Jessa, adapted from Claudia Kleinman
#'
#' @examples
#' mm2hg("Tubb5", "Nes")
#' mm2hg(c("Tubb5", "Nes"))
#'
#' mat1 <- matrix(seq(1, 15), nrow = 3)
#' rownames(mat1) <- c("Gfap", "Fabp7", "Foo")
#' mm2hg(mat1)
mm2hg <- function (x, ...) {
    UseMethod("mm2hg", x)
}

#' hg2mm
#'
#' Convert genes from human to mouse, wrapper around \code{\link{convertGenesBetweenSpeciesHomologs}}.
#' For vectors, the output is \strong{not order preserving}.
#' For matrices, order is preserved, but \strong{genes where the homolog is not
#' found are silently dropped}.
#'
#' @param x Character vector or matrix with human gene symbols to convert.
#' If character, provide either individual strings or one character vector.
#' If a matrix, rownames should be gene symbols.
#'
#' @export
#' @author Selin Jessa, adapted from Claudia Kleinman
#'
#' @examples
#' hg2mm("TUBB", "NES")
#' hg2mm(c("TUBB", "NES"))
hg2mm <- function (x, ...) {
    UseMethod("hg2mm", x)
}


#' @export
convertGenesBetweenSpeciesHomologs.character <- function(symbol, from, to = NULL) {
    switch(from, "hg" = hg2mm.genes, "mm" = mm2hg.genes) %>%
        filter(toupper(gene_symbol) %in% toupper(symbol)) %>%
        .$homologous_gene_symbol
}


#' @export
mm2hg.character <- function(...) {
    convertGenesBetweenSpeciesHomologs(c(...), "mm")
}

#' @export
hg2mm.character <- function(...) {
    convertGenesBetweenSpeciesHomologs(c(...), "hg")
}


#' @export
convertGenesBetweenSpeciesHomologs.matrix <- function(x, from, to = NULL) {

    annotation <- switch(from, "hg" = hg2mm.genes, "mm" = mm2hg.genes)

    # 1. Make a gene map from one species to the other
    gene_map <- annotation %>%
        # TODO: This drops gene names if there are duplicates
        distinct(gene_symbol, .keep_all = TRUE) %>%
        filter(!is.na(gene_symbol), !is.na(homologous_gene_symbol)) %>%
        select(gene_symbol, homologous_gene_symbol) %>%
        as.data.frame() %>%
        tibble::column_to_rownames(., var = "gene_symbol")

    # 2. Convert the rownames on the count matrix
    converted <- gene_map[rownames(x), ]
    found <- which(!is.na(converted))
    x_subset <- x[found, ]
    rownames(x_subset) <- gene_map[rownames(x_subset), ]

    return(x_subset)

}

#' @export
mm2hg.matrix <- function(x, ...) {

    convertGenesBetweenSpeciesHomologs.matrix(x = x, from = "mm", to = "hg")

}


#' @export
hg2mm.matrix <- function(x, ...) {

    convertGenesBetweenSpeciesHomologs.matrix(x = x, from = "hg", to = "mm")

}



#' filterGenesSample
#'
#' Given a set of genes, filter to the ones found in the sample.
#' This is useful because some Seurat functions (Seurat::VlnPlot, Seurat::FeaturePlot)
#' will quit on an error if a gene is not found (i.e. not in the rownames of seurat@@data).
#'
#' @param genes Character vector of genes to search for
#' @param seurat Seurat object
#' @param from_sp One of "hg" or "mm". Default: "hg"
#' @param to_sp One of "hg" or "mm". Default: "mm"
#'
#' @return Character vector of converted gene names which are detected in the sample
#' @export
#' @author adapted from Claudia Kleinman
#'
#' @examples
#' filterGenesSample(c("IL32", "foo"), pbmc, "hg", "hg")
filterGenesSample <- function(genes, seurat, from_sp = "hg", to_sp = "mm") {

    sample_genes <- rownames(seurat@hvg.info)
    if (from_sp == to_sp) {
        filter(geneName(genes, from_sp), gene_symbol %in% sample_genes) %>% .$gene_symbol
    }
    else {
        conv_genes <- convertGenesBetweenSpeciesHomologs(genes, from_sp)
        conv_genes[which(conv_genes %in% sample_genes)]
    }

}


#' Convert between gene symbols and ENSEMBL gene IDs
#'
#' @param x Character vector of genes, or matrix (e.g. of counts)
#' where rownames are genes.
#' @param sp Character, one of "mm" or "hg", specifying the species for which
#' data is given in \code{mat}. Default: "hg".
#' @param from Character vector, one of "gene_symbol" or "ensembl_gene_id",
#' specifying in which format the genes in \code{x} are provided.
#'
#' @author Selin Jessa
convertGenesBetweenSymbolsAndID <- function(x, sp = "hg", from = "gene_symbol") {

    annotation <- switch(sp, "hg" = hg.genes, "mm" = mm.genes)

    # 1. Make a map from gene name to ENSEMBL id
    gene_map <- annotation %>%
        # TODO: This drops gene names if there are duplicates
        distinct(gene_symbol, .keep_all = TRUE) %>%
        # Filter out things that we can't map
        filter(!is.na(gene_symbol)) %>%
        filter(!is.na(ensembl_gene_id)) %>%
        select(gene_symbol, ensembl_gene_id) %>%
        as.data.frame() %>%
        tibble::column_to_rownames(var = from)

    # 2. Convert the rownames on the count matrix
    if (is.character(x)) {

        converted <- gene_map[x, ]
        found <- which(!is.na(converted))
        x_found <- x[found]
        x_converted <- gene_map[x_found, ]

    } else if (is.matrix(x)) {

        converted <- gene_map[rownames(x), ]
        found <- which(!is.na(converted))
        x_converted <- x[found, ]
        rownames(x_converted) <- gene_map[rownames(x_converted), ]

    }

    return(x_converted)

}

#' Convert ENSEMBL gene IDs to gene symbols
#'
#' @param x Character vector of genes, or matrix (e.g. of counts)
#' where rownames are genes.
#' @param sp Character, one of "mm" or "hg", specifying the species for which
#' data is given in \code{mat}. Default: "hg".
#'
#' @export
#' @author Selin Jessa
#'
#' @examples
#' mat1 <- matrix(seq(1, 15), nrow = 3)
#' rownames(mat1) <- c("ENSMUSG00000020932", "ENSMUSG00000019874", "Foo")
#' ensembl2symbols(mat1)
#'
#' ensembl2symbols(c("ENSMUSG00000020932", "ENSMUSG00000019874"), "mm")
ensembl2symbols <- function(x, sp = "hg") {

    convertGenesBetweenSymbolsAndID(x = x, sp = sp, from = "ensembl_gene_id")

}

#' Convert gene symbols to ENSEMBL gene IDs
#'
#' @param x Character vector of genes, or matrix (e.g. of counts)
#' where rownames are genes.
#' @param sp Character, one of "mm" or "hg", specifying the species for which
#' data is given in \code{mat}. Default: "hg".
#'
#' @export
#' @author Selin Jessa
#'
#' @examples
#' mat1 <- matrix(seq(1, 15), nrow = 3)
#' rownames(mat1) <- c("Gfap", "Fabp7", "Foo")
#' symbols2ensembl(mat1)
#'
#' symbols2ensembl(c("Gfap", "Fabp7"), "mm")
symbols2ensembl <- function(x, sp = "hg") {

    convertGenesBetweenSymbolsAndID(x = x, sp = sp, from = "gene_symbol")

}


#' Annotate cluster markers computed by Seurat
#'
#' @param markers A dataframe containing cluster markers as computed by
#' \link[Seurat]{FindAllMarkers}. Rownames should be gene symbols.
#' @param species String specifying the species, one of "mm" or "hg", used to determine the
#' annotation. Default: "mm".
#'
#'
#' @author Selin Jessa
#' @export
#' @examples
#' annotateMarkers(markers_pbmc, "hg")
annotateMarkers <- function(markers, species = "mm") {

    annotation <- switch(species, "hg" = hg.genes, "mm" = mm.genes)

    markers %>%
        tibble::rownames_to_column(var = "external_gene_name") %>%
        left_join(annotation, by = c("external_gene_name" = "gene_symbol")) %>%
        distinct(cluster, external_gene_name, .keep_all = TRUE) %>%
        select(cluster, external_gene_name, avg_logFC, p_val_adj, pct.1, pct.2, description, p_val, gene_biotype, ensembl_gene_id) %>%
        arrange(cluster, desc(avg_logFC))


}

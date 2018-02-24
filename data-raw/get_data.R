# Selin Jessa

# Get data from Seurat to use for testing and examples

pbmc <- Seurat::pbmc_small
markers_pbmc <- Seurat::FindAllMarkers(pbmc)

usethis::use_data(pbmc)
usethis::use_data(markers_pbmc)

# Annotations
# Obtained these annotations from Claudia

hg.genes <- fread("~/scratch/cortex/annotations/hg_gene_annotation.txt", data.table = FALSE)
mm.genes <- fread("~/scratch/cortex/annotations/mm_gene_annotation.txt", data.table = FALSE)
mm2hg.genes <- fread("~/scratch/cortex/annotations/mm_gene_annotation_hg_homol.txt", data.table = FALSE)
hg2mm.genes <- fread("~/scratch/cortex/annotations/hg_gene_annotation_mm_homol.txt", data.table = FALSE)

usethis::use_data(hg.genes)
usethis::use_data(mm.genes)
usethis::use_data(mm2hg.genes)
usethis::use_data(hg2mm.genes)

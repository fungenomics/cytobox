# Get data from Seurat to use for testing and examples

pbmc <- Seurat::pbmc_small
markers_pbmc <- Seurat::FindAllMarkers(pbmc)

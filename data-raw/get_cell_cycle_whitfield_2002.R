# Alexis Blanchet-Cohen

# Data obtained from Whitfield et al, 2012
library(data.table)
library(dplyr)
library(usethis)

#cell.cycle.genes.whitfield.2002 <- fread("/gs/project/txk-653-aa/blancha/analyses/misc/CellCycleGeneList_1134_whitfield_2002.txt", data.table = FALSE)
#cell.cycle.genes.whitfield.2002 <- select(cell.cycle.genes.whitfield.2002, -CLONEID)

#cell.cycle.genes.whitfield.2002$PHASE <- gsub("*", "", cell.cycle.genes.whitfield.2002$PHASE, fixed=TRUE)
#cell.cycle.genes.whitfield.2002 <- arrange(cell.cycle.genes.whitfield.2002, PHASE)
#cell.cycle.genes.whitfield.2002$gene_symbol <- sapply(strsplit(cell.cycle.genes.whitfield.2002$NAME, split=" ", fixed= TRUE), "[", 1)

#cell.cycle.genes.whitfield.2002 <- select(cell.cycle.genes.whitfield.2002, gene_symbol, PHASE)
#colnames(cell.cycle.genes.whitfield.2002) <- tolower(colnames((cell.cycle.genes.whitfield.2002)))

cell.cycle.genes.whitfield.2012 <- fread("/gs/project/txk-653-aa/blancha/analyses/njabado/single_cell/samples/ET_CT_12/data/CellCycleGeneList_1134_whitfield_2002_mice_gene_symbols.txt", data.table=FALSE)

colnames(cell.cycle.genes.whitfield.2012) <- tolower(colnames(cell.cycle.genes.whitfield.2012))
colnames(cell.cycle.genes.whitfield.2012)[1] <- "mmusculus.gene.symbol"

usethis::use_data(cell.cycle.genes.whitfield.2012, overwrite=TRUE)

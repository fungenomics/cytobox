# Alexis Blanchet-Cohen

# Abcam neuronal markers annotations obtained from http://www.abcam.com/neuroscience/neural-markers-guide.
# Manually formatted and annotated to correspond to Matt's database format.
library(data.table)
library(usethis)

abcam.neuronal.markers <- fread("../../../analyses/njabado/single_cell/abcam_neural_markers/abcam_neuronal_markers.txt", data.table = FALSE)

usethis::use_data(abcam.neuronal.markers, overwrite=TRUE)

install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)

seu <- readRDS("data/processed/AMC_subset_annotated.rds")
SaveH5Seurat(seu, filename = "data/processed/AMC_subset_annotated.h5Seurat")
Convert("data/processed/AMC_subset_annotated.h5Seurat", dest = "h5ad")

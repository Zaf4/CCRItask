library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
setwd("~/dev/CCRItask")

# List all H5 files in the raw data directory
h5_files <- list.files("data/raw/GSE147821_RAW", pattern = ".h5$", full.names = TRUE)

# map GSM ids to info
files <- list(
  "GSM4446535" = "week8_001",
  "GSM4446536" = "week9_063",
  "GSM4446537" = "week6_088",
  "GSM4446538" = "week14_123",
  "GSM4446539" = "week12_124",
  "GSM4446540" = "week8_125",
  "GSM4446541" = "week9_005",
  "GSM4446542" = "week11_006",
  "GSM4446543" = "week9_007",
  "GSM4734601" = "week8_016",
  "GSM4734602" = "week9_031_paraganglia",
  "GSM4734603" = "week12_035",
  "GSM4734604" = "week12_036_extraadrenal"
)

# Load each sample into a list of Seurat objects
seurat_list <- lapply(h5_files, function(file) {
  sample_name <- gsub(".*GSM\\d+_10X_\\d+_(\\d+).*", "\\1", basename(file))
  print(sample_name)
  seurat_obj <- CreateSeuratObject(
    counts = Read10X_h5(file),
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  seurat_obj$sample <- sample_name
  return(seurat_obj)
})


seurat_list <- lapply(names(h5_files), function(gsm_id) {
  # Get the sample name from our files dictionary
  sample_name <- files[[gsm_id]]
  
  # Extract week number (digits after "week")
  week <- gsub("week(\\d+).*", "\\1", sample_name)
  
  # Extract sample ID (digits after underscore)
  sample_id <- gsub(".*_(\\d+).*", "\\1", sample_name)
  
  # Read the data and create Seurat object
  seurat_obj <- CreateSeuratObject(
    counts = Read10X_h5(h5_files[[gsm_id]]),
    project = sample_name,
    min.cells = 3,
    min.features = 200
  )
  
  # Add metadata
  seurat_obj$sample <- sample_name
  seurat_obj$week <- week
  seurat_obj$sample_id <- sample_id
  seurat_obj$gsm_id <- gsm_id  # Store the GSM ID as well
  
  return(seurat_obj)
})


# Name the list elements (optional)
names(seurat_list) <- sapply(seurat_list, function(x) unique(x$sample))


# QC: Filter low-quality cells
seurat_list <- lapply(seurat_list, function(x) {
  x[["percent.mt"]] <- PercentageFeatureSet(x, pattern = "^MT-")
  x <- subset(x, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 15)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# Add cell cycle correction
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
seurat_list <- lapply(seurat_list, function(x) {
  x <- CellCycleScoring(
    x,
    s.features = s_genes,
    g2m.features = g2m_genes,
    set.ident = FALSE
  )
  x$CC.Difference <- x$S.Score - x$G2M.Score  # Useful for regression
  return(x)
})

# scaling and regression for cell cycle
seurat_list <- lapply(seurat_list, function(x) {
  x <- ScaleData(
    x,
    vars.to.regress = "CC.Difference",  
    verbose = FALSE
  )
  return(x)
})

# Select integration anchors
anchors <- FindIntegrationAnchors(
  object.list = seurat_list,
  dims = 1:30,
  anchor.features = 2000,
  normalization.method = "LogNormalize"
)

# Integrate data
integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:30,
  new.assay.name = "integrated"
)

# Switch to integrated assay for downstream analysis
DefaultAssay(integrated) <- "integrated"

# Scale data and run PCA
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)

# Check batch effects (should be mixed in PCA)
DimPlot(integrated, reduction = "pca", group.by = "sample")

# Run UMAP and cluster
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

# Visualize clusters (colored by cluster or sample)
DimPlot(integrated, group.by = "seurat_clusters", label = TRUE)
DimPlot(integrated, group.by = "sample")
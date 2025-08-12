rm(list=ls())
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(harmony)

samples = read_delim('Penter_data/download_samples.csv') %>% 
  pull(title)

x = lapply(samples, function(s) {
  dd = Read10X(paste0('Penter_data/', s))  
  # CreateSeuratObject(counts = dd, project = s, min.cells = 3, min.features = 200)
  return(dd)
})
names(x) = samples

# create seurat objects
xs = lapply(x, function(s) {CreateSeuratObject(s)})
xs = lapply(xs, function(s) {
  # The [[ operator can add columns to object metadata. This is a great place to stash QC stats
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  return(s)
})
saveRDS(xs, 'Hemasphere_data/penter_seurat.rds')

# Visualize QC metrics as a violin plot
VlnPlot(xs$CLL1_3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

################################################################################################

# do some filtering and QC
# Thresholds
min_lib_size <- 900
min_n_genes <- 250
max_pct_mt <- 22.5
min_cells <- 3

seurat_list = lapply(xs %>% names, function(x) {
  
  print(x)
  
  x = xs[[x]]
  
  metadata_before_qc = x@meta.data
  
  is_low_quality <- 
    x$nCount_RNA < min_lib_size |
    x$nFeature_RNA < min_n_genes |
    x$percent.mt > max_pct_mt
  
  x <- subset(x, cells = colnames(x)[!is_low_quality])
  
  kept_genes <- rownames(x)[n_cells > min_cells]
  x = subset(x, features = kept_genes)
  
  return(x)
})
names(seurat_list) = names(xs)
saveRDS(seurat_list, 'Hemasphere_data/penter_seurat_filtered_qc.rds')
################################################################################################

color_palette <- c("black", "gray", "red", "yellow", "violet", "green4",
                   "blue", "mediumorchid2", "coral2", "blueviolet",
                   "indianred4", "deepskyblue1", "dimgray", "deeppink1",
                   "green", "lightgray", "hotpink1")

seurat_list = seurat_list[-7]
seurat_list_v2 = lapply(seurat_list %>% names, function(s) {
  
  print(s)
  seurat = seurat_list[[s]]
  seurat$sample_id = rep(s, seurat$orig.ident %>% length)
  
  seurat <- seurat %>%
    NormalizeData(
      normalization.method = "LogNormalize",
      scale.factor = 10000
    ) %>%
    FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>% 
    ScaleData(features = rownames(seurat)) %>%
    RunPCA() %>%
    # RunHarmony(group.by.vars = "sample_id", # reduction = "pca", 
    #            # dims = 1:30
    #            ) #%>%
    RunUMAP(dims = 1:30, reduction = "pca")
  return(seurat)
})
names(seurat_list_v2) = names(seurat_list)

saveRDS(seurat_list_v2, 'Hemasphere_data/penter_seurat_clustered.rds')

# lineage_markers <- c("CD79A", "CD3D",  "NKG7",
#                      "LYZ", "TOP2A")
# lineage_markers_ggs <- purrr::map(lineage_markers, function(x) {
#   p <- FeaturePlot(seurat_list_v2$CLL5_1, features = x, pt.size = 0.45)
#   p +
#     scale_color_viridis_c(option = "magma")
# })
# names(lineage_markers_ggs) <- lineage_markers
# lineage_markers_ggs

# (FeaturePlot(seurat_list_v2$CLL1_3.1, "nFeature_RNA") +
#   scale_color_viridis_c(option = "magma")) + 
# (FeaturePlot(seurat_list_v2$CLL1_3.1, "percent.mt") +
#      scale_color_viridis_c(option = "magma")) 

seurat_clustered = lapply(seurat_list_v2, function(seurat) {
  seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:20)
  seurat <- FindClusters(seurat, resolution = 0.3)
  return(seurat)
})

wrap_plots(lapply(
  seurat_clustered, DimPlot
))

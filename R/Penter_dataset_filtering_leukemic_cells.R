rm(list=ls())
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(patchwork)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/TLS/CLL_revision')

seurat_list = readRDS('Hemasphere_data/penter_seurat_filtered_qc.rds')
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

lineage_markers <- c("CD79A", "CD3D", "NKG7", "LYZ", "TOP2A")
lineage_markers_ggs <- lapply(seurat_list_v2, function(s) {
  pp = purrr::map(lineage_markers, function(x) {
    p <- FeaturePlot(s, features = x, pt.size = 0.45)
    p +
      scale_color_viridis_c(option = "magma")
  })
  names(pp) <- lineage_markers
  return(pp)
})

pdf('Hemasphere_data/Penter_samples_lineage_markers.pdf', width = 12, height = 8)  
lapply(lineage_markers_ggs %>% names, function(x) {
  wrap_plots(lineage_markers_ggs[[x]]) + 
    plot_annotation(title = x)
})
dev.off()


# (FeaturePlot(seurat_list_v2$CLL1_3.1, "nFeature_RNA") +
#   scale_color_viridis_c(option = "magma")) + 
# (FeaturePlot(seurat_list_v2$CLL1_3.1, "percent.mt") +
#      scale_color_viridis_c(option = "magma")) 

seurat_clustered = lapply(seurat_list_v2, function(seurat) {
  seurat <- FindNeighbors(seurat, reduction = "pca", dims = 1:20)
  seurat <- FindClusters(seurat, resolution = 0.1)
  return(seurat)
})

pdf('Hemasphere_data/penter_clustering.pdf')
lapply(
  seurat_clustered %>% names, function(x) {
    DimPlot(seurat_clustered[[x]], cols = color_palette) + 
      ggtitle(x)
  }
)
dev.off()

# trying to find subcluster in specific groups
seurat_clustered$CLL2_1 = FindSubCluster(
  seurat_clustered$CLL2_1,
  cluster = "2",
  graph.name = "RNA_snn",
  subcluster.name = "fetch_T_cells",
  resolution = 0.2
)
DimPlot(seurat_clustered$CLL2_1, cols = color_palette, group.by = 'fetch_T_cells')

seurat_clustered$CLL3_1 = FindSubCluster(
  seurat_clustered$CLL3_1,
  cluster = "1",
  graph.name = "RNA_snn",
  subcluster.name = "fetch_T_cells",
  resolution = 0.2
)
DimPlot(seurat_clustered$CLL3_1, cols = color_palette, group.by = 'fetch_T_cells')

seurat_clustered$CLL4_1 = FindSubCluster(
  seurat_clustered$CLL4_1,
  cluster = c("1", '2'),
  graph.name = "RNA_snn",
  subcluster.name = "fetch_T_cells",
  resolution = 0.2
)
DimPlot(seurat_clustered$CLL4_1, cols = color_palette, group.by = 'fetch_T_cells')

seurat_clustered$CLL5_1 = FindSubCluster(
  seurat_clustered$CLL5_1,
  cluster = c("1", '2'),
  graph.name = "RNA_snn",
  subcluster.name = "fetch_T_cells",
  resolution = 0.2
)
DimPlot(seurat_clustered$CLL5_1, cols = color_palette, group.by = 'fetch_T_cells')

seurat_clustered$CLL6_1 = FindSubCluster(
  seurat_clustered$CLL6_1,
  cluster = '1',
  graph.name = "RNA_snn",
  subcluster.name = "fetch_T_cells",
  resolution = 0.2
)
DimPlot(seurat_clustered$CLL6_1, cols = color_palette, group.by = 'fetch_T_cells')

seurat_clustered$CLL8_1 = FindSubCluster(
  seurat_clustered$CLL8_1,
  cluster = c('1', '3'),
  graph.name = "RNA_snn",
  subcluster.name = "fetch_T_cells",
  resolution = 0.2
)
DimPlot(seurat_clustered$CLL8_1, cols = color_palette, group.by = 'fetch_T_cells')


pdf('Hemasphere_data/penter_clustering_t_cells.pdf')
lapply(
  seurat_clustered[-1] %>% names, function(x) {
    DimPlot(seurat_clustered[[x]], cols = color_palette, group.by = 'fetch_T_cells') + 
      ggtitle(x)
  }
)
dev.off()

# leukemic cells
# CLL1_3.1	C0,C1,C2
# CLL2_1	C0,C1
# CLL3_1	C0, C1_2, C2
# CLL4_1	C0
# CLL5_1	C0,C3, C2_0, C1_3
# CLL6_1	C0,C3, C1_2
# CLL8_1	C0, C2

# select only leukemic clusters!
leuk_clusters = list(
  'CLL1_3.1' = c('0','1','2'), 
  'CLL2_1' = c('0','1'), 
  'CLL3_1' = c('0', '1_2', '2'), 
  'CLL4_1' = c('0'), 
  'CLL5_1' = c('0', '1_3', '2_0', '3'), 
  'CLL6_1' = c('0', '1_2', '3'), 
  'CLL8_1' = c('0', '2')
)

leukemic_cells = lapply(names(leuk_clusters)[-1], function(x) {
  Idents(seurat_clustered[[x]]) = 'fetch_T_cells'
  cl = leuk_clusters[[x]]
  subset(seurat_clustered[[x]], idents = cl)
}) 
names(leukemic_cells) = names(leuk_clusters)[-1]
leukemic_cells = c(
  list('CLL1_3.1' = subset(seurat_clustered$CLL1_3.1, idents = leuk_clusters$CLL1_3.1)), 
  leukemic_cells
)
saveRDS(leukemic_cells, 'Hemasphere_data/penter_seurat_leukemic_clusters.rds')

dir.create('Hemasphere_data/leukemic_penter')
lapply(names(leukemic_cells), function(x) {
  saveRDS(leukemic_cells[[x]], paste0('Hemasphere_data/leukemic_penter/', x, '.rds'))
})



rm(list=ls())
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(patchwork)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/TLS/CLL_revision')

setwd('Hemasphere_data/leukemic_penter')
leukemic_cells = lapply(list.files(), readRDS)
names(leukemic_cells) = gsub('.rds', '', list.files())

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/TLS/CLL_revision')

process_seurat <- function(seurat_obj, dims = 1:30) {
  
  # Idents(seurat_obj) = 'fetch_T_cells'
  
  seurat_obj %>%
    FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>% 
    ScaleData(features = rownames(seurat_obj)) %>%
    RunPCA() %>%
    RunUMAP(dims = dims, reduction = "pca")
}

seurat_list = lapply(leukemic_cells, function(x){
  
  if(x$sample_id %>% unique == 'CLL1_3.1') {
    Idents(x) = 'seurat_clusters'
  } else {
    Idents(x) = 'fetch_T_cells'
  }
  
  process_seurat(x, dims = 1:20)
})

seurat_list = lapply(seurat_list, function(s) {
  Idents(s) = 'fetch_T_cells'
  return(s)
})

Idents(seurat_list$CLL1_3.1) = 'seurat_clusters'

pdf('Hemasphere_data/penter_leukemic_clusters_expr_v2.pdf')
lapply(seurat_list %>% names, function(s) {
  mt = FeaturePlot(seurat_list[[s]], "percent.mt") +
    scale_color_viridis_c(option = "magma")
  rna = FeaturePlot(seurat_list[[s]], "nFeature_RNA") +
    scale_color_viridis_c(option = "magma")
  
  (mt + rna) + 
    plot_annotation(title = s)
})
dev.off()

seurat_list = lapply(seurat_list, function(x) {
  # Idents(x) = 'fetch_T_cells'
  x <- FindNeighbors(x, reduction = "pca", dims = 1:20)
  x <- FindClusters(x, resolution = 0.3)
  return(x)
})

pdf('Hemasphere_data/leukemic_new_clusters_v2.pdf')
lapply(seurat_list, DimPlot)
dev.off()

# lapply(seurat_list, function(s) {
#   Idents(s) = 'fetch_T_cells'
# })
# see if there are any clusters to remove from the objects

# remove clusters with high mytocondrial expression and low mumber of features
seurat_list$CLL1_3.1 = subset(seurat_list$CLL1_3.1, seurat_clusters != '2')
seurat_list$CLL3_1 = subset(seurat_list$CLL3_1, seurat_clusters != '5')
seurat_list$CLL8_1 = subset(seurat_list$CLL8_1, seurat_clusters != '1')

seurat_list$CLL1_3.1 = process_seurat(seurat_list$CLL1_3.1, dims = 1:20)
seurat_list$CLL3_1 = process_seurat(seurat_list$CLL3_1, dims = 1:20)
seurat_list$CLL8_1 = process_seurat(seurat_list$CLL8_1, dims = 1:20)

seurat_list$CLL1_3.1 = seurat_list$CLL1_3.1 %>% 
  FindNeighbors(., reduction = 'pca', dims = 1:20) %>% 
  FindClusters(., resolution = 0.3)
seurat_list$CLL3_1 = seurat_list$CLL3_1 %>% 
  FindNeighbors(., reduction = 'pca', dims = 1:20) %>% 
  FindClusters(., resolution = 0.3)
seurat_list$CLL8_1 = seurat_list$CLL8_1 %>% 
  FindNeighbors(., reduction = 'pca', dims = 1:20) %>% 
  FindClusters(., resolution = 0.3)

pdf('Hemasphere_data/CCND2_expression_penter.pdf')
lapply(seurat_list, function(x) {
  FeaturePlot(x, "CCND2") + scale_color_viridis_c(option = "magma")
})
dev.off()

# looking for other microenvironment contaminations
microenv_markers <- c("CD3D", "NKG7", "LYZ")

pdf('Hemasphere_data/residual_microenvironment_contamination.pdf')

lapply(seurat_list %>% names, function(s) {
  
  ss = seurat_list[[s]]
  
  feat_plots_microenv <- purrr::map(microenv_markers, function(x) {
    
    p <- FeaturePlot(ss, x) +
      scale_color_viridis_c(option = "magma") 
    p
  })
  
  wrap_plots(feat_plots_microenv) + 
    plot_annotation(title = s)
  
})
dev.off()

# now annotate cells!
color_palette <- c("black", "gray", "red", "yellow", "violet", "green4",
                   "blue", "mediumorchid2", "coral2", "blueviolet",
                   "indianred4", "deepskyblue1", "dimgray", "deeppink1",
                   "green", "lightgray", "hotpink1")


# compute cell cycle scoring, since it is an indication of RT
seurat_list = lapply(seurat_list, function(seurat) {
  CellCycleScoring(
    seurat,
    s.features = cc.genes.updated.2019$s.genes,
    g2m.features = cc.genes.updated.2019$g2m.genes,
    set.ident = FALSE
  )
})

DimPlot(seurat_list$CLL3_1, group.by = 'Phase')

# get markers to annotate the cell clusters
min_log2FC <- 0.3
alpha <- 0.001

markers = lapply(seurat_list, function(seurat) {
  m <- FindAllMarkers(seurat, only.pos = TRUE, logfc.threshold = min_log2FC)
  m %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(p_val_adj < alpha) %>%
    arrange(cluster) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC), .by_group = TRUE)
})
saveRDS(markers, 'Hemasphere_data/markers_results_new.rds')

# CXCR4 	CXCR4hiCD27lo
# CD27 	CXCR4loCD27hi
# CCND2 	CCND2hi RT
# PAGE2, PAGE2B 	CCND2lo RT
# TOP2A, MKI67 	RT proliferative
# MIR155HG 	MIR155HGhi
# IGHM, MZB1, XBP1 	MZB1hiIGHMhiXBP1hi
markers = readRDS('Hemasphere_data/markers_results_new.rds')
mm = lapply(markers, function(m) {
  m %>%
    dplyr::filter(gene %in% c('CXCR4', 'CD27', 'CCND2', 'PAGE2', 'PAGE2B', 'TOP2A', 'MKI67', 'MIR155HG', 'IGHM', 'MZB1', 'XBP1'))
})

# DimPlot(seurat_list$CLL3_1)
# FeaturePlot(seurat_list$CLL3_1, features = 'CXCR4')

# cll1.3
# seurat_list$CLL1_3.1$annotation_final <- factor(
#   seurat_list$CLL1_3.1$RNA_snn_res.0.3,
#   levels = c("0", "3", "2", '3'),
# )
# new_levels_12 <- c("CXCR4hiCD27lo", "CXCR4loCD27hi", "other_1", 'other2')
# levels(seurat_list$CLL1_3.1$annotation_final) <- new_levels_12
# 
# Idents(seurat_list$CLL1_3.1) <- "annotation_final"

# cll2
seurat_list$CLL2_1$annotation_final <- factor(
  seurat_list$CLL2_1$RNA_snn_res.0.3,
  levels = c("0", "2", "1", '3'),
)
new_levels_12 <- c("CXCR4hiCD27lo", "CXCR4loCD27hi", "other", 'other2')
levels(seurat_list$CLL2_1$annotation_final) <- new_levels_12

Idents(seurat_list$CLL2_1) <- "annotation_final"

# cll3
seurat_list$CLL3_1$annotation_final <- factor(
  seurat_list$CLL3_1$RNA_snn_res.0.3,
  levels = c("3", "1", "2", '0'),
)
new_levels_12 <- c("CXCR4hiCD27lo", "CXCR4loCD27hi", "other1", 'other2')
levels(seurat_list$CLL3_1$annotation_final) <- new_levels_12

Idents(seurat_list$CLL3_1) <- "annotation_final"

# cll6
seurat_list$CLL6_1$annotation_final <- factor(
  seurat_list$CLL6_1$RNA_snn_res.0.3,
  levels = c("2", "1", "3", '4', '5', '0'),
)
new_levels_12 <- c("CXCR4hiCD27lo", "CXCR4loCD27hi", "other1", 'other2', 'other3', 'other4')
levels(seurat_list$CLL6_1$annotation_final) <- new_levels_12

Idents(seurat_list$CLL6_1) <- "annotation_final"

good_samples = seurat_list[c("CLL2_1", "CLL3_1", "CLL6_1")]
# saveRDS(good_samples, 'Hemasphere_data/good_samples_clustered_annotated.rds')

# ppp = c(
#   'ALDOA',	
#   'ALDOB',	
#   'ALDOC',	
#   'DERA',	
#   'FBP1',	
#   # 'FBP2',	
#   'G6PD',	
#   'GPI',	
#   'H6PD',	
#   'PFKL',	
#   'PFKM',	
#   'PFKP',	
#   'PGD',	
#   'PGLS',	
#   'PGM1',	
#   'PGM2',	
#   'PRPS1',	
#   # 'PRPS1L1',	
#   'PRPS2',	
#   'RBKS',	
#   'RPE',	
#   # 'RPEL1',	
#   'RPIA',	
#   'TALDO1',	
#   'TKT',	
#   'TKTL1',	
#   'TKTL2'
# )

ppp = c(
  'G6PD', 
  'PGLS', 
  'PGD', 
  'RPE', 
  'RPIA', 
  'TKT', 
  'TALDO1'
)

# create dotplots 


lapply(good_samples %>% names, function(x) {
  
  pdf(paste0('Hemasphere_data/ppp_different_populations_', x, '.pdf'), width = 9, height = 5)
  print(DotPlot(object = good_samples[[x]], features = c(ppp, 'CXCR4', 'CD27', 'CD5'), cols = 'RdBu',
          idents = c('CXCR4hiCD27lo', 'CXCR4loCD27hi'), dot.scale = 19) +
    theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1, size = 16), 
          legend.position = 'right', legend.direction = 'vertical', 
          axis.line = element_line(linewidth = 0.3), 
          axis.ticks = element_line(linewidth = 0.3)) +
    xlab('') + 
    ylab('')) 
  dev.off()
})

saveRDS(good_samples, 'Hemasphere_data/good_samples_penter.rds')

lapply(good_samples, function(s){
  
  cells_to_vis = s@meta.data %>% 
    dplyr::filter(annotation_final %in% c('CXCR4hiCD27lo', 'CXCR4loCD27hi')) %>% 
    tibble::rownames_to_column('id') %>% 
    group_by(annotation_final) %>% 
    group_split() 
  
  names(cells_to_vis) = lapply(cells_to_vis, function(x) {
    x %>% 
      pull(annotation_final) %>% unique
  }) %>% unlist
  
  cells_to_vis = lapply(cells_to_vis, function(x) {x$id})
  
  sid = s$sample_id %>% unique
  
  pdf(paste0('Hemasphere_data/', sid, '_umap.pdf'))
  print(DimPlot(s, cells.highlight = cells_to_vis, cols.highlight = c('CXCR4hiCD27lo' = '#16325B', 'CXCR4loCD27hi' = '#0F828C')) + 
    theme(legend.position = 'bottom') + 
    ggtitle(sid))
  dev.off()
})


# FeaturePlot(object = good_samples$, features = 'G6PD')

# perform differential expression analysis
run_dge = function(x, idents, name_ident) {
  Idents(x) <- name_ident
  dge = FindMarkers(x, 
                    ident.1 = idents[1], 
                    ident.2 = idents[2], 
                    only.pos = FALSE,
                    logfc.threshold = 0)
  return(dge)
}

dge = lapply(good_samples, function(x){
  # lapply(x, function(s) {
    run_dge(x, idents = rev(c('CXCR4hiCD27lo', 'CXCR4loCD27hi')), name_ident = 'annotation_final')
  # })
})
saveRDS(dge, 'Hemasphere_data/dge_penter_all.rds')

ppp_dge_penter = lapply(dge %>% names, function(x) {
  dge[[x]] %>% 
    dplyr::mutate(sample = x) %>% 
    tibble::rownames_to_column('gene') %>% 
    dplyr::filter(gene %in% c(ppp, 'CXCR4', 'CD27', 'CD5'))
}) %>% bind_rows()

write.table(ppp_dge_penter, 'Hemasphere_data/ppp_dge_penter.csv', sep = ',', quote = F, row.names = F, col.names = T)

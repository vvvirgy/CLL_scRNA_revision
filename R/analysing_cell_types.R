rm(list=ls())

.libPaths()
library(tidyverse)
library(Seurat)
library(SingleCellComplexHeatMap)

data = '/orfeo/LTS/CDSLab/LT_storage/vgazziero/cll_tls_prj/Nadeu2022_NatMed_scRNAseq_data'
list.files(data)
seurat_obj = paste0(data, '/seurat_objects')
list.files(seurat_obj)

objs = list(
  's19' = readRDS(paste0(seurat_obj, '/6.seurat_annotated_19.rds')), 
  's3299' = readRDS(paste0(seurat_obj, '/6.seurat_annotated_3299.rds'))
)

objs_by_time = lapply(objs, function(x) {
  Seurat::SplitObject(x, split.by = 'time_point')
})


# DimPlot(objs$s19, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = 'time_point') 
# DimPlot(objs$s3299, reduction = "umap", label = TRUE, pt.size = 0.5, group.by = "time_point") 

ppp = c(
 'ALDOA',	
 'ALDOB',	
 'ALDOC',	
 'DERA',	
 'FBP1',	
 'FBP2',	
 'G6PD',	
 'GPI',	
 'H6PD',	
 'PFKL',	
 'PFKM',	
 'PFKP',	
 'PGD',	
 'PGLS',	
 'PGM1',	
 'PGM2',	
 'PRPS1',	
 'PRPS1L1',	
 'PRPS2',	
 'RBKS',	
 'RPE',	
 'RPEL1',	
 'RPIA',	
 'TALDO1',	
 'TKT',	
 'TKTL1',	
 'TKTL2'
)

pdf('Hemasphere_data/s19_ppp.pdf', width = 20, height = 20)
lapply(objs_by_time$s19 %>% names, function(x) {
  VlnPlot(objs_by_time$s19[[x]], features = ppp) + 
    patchwork::plot_annotation(title = x)
})
dev.off()

pdf('Hemasphere_data/s3299_ppp.pdf', width = 20, height = 20)
lapply(objs_by_time$s3299 %>% names, function(x) {
  VlnPlot(objs_by_time$s3299[[x]], features = ppp) + 
    patchwork::plot_annotation(title = x)
})
dev.off()

DotPlot(object = objs_by_time$s19$T3, features = c('TKT', 'CXCR4'), 
        # idents = c('CXCR4hiCD27lo', 'CXCR4loCD27hi', 'RT proliferative')
        ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
ggsave('Hemasphere_data/tkt_expression.png', bg = 'white')

DotPlot(object = objs_by_time$s3299$T1, features = 'CXCR4', cols = 'YlGn') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))





df = LayerData(objs_by_time$s19$T3[c('TKT', 'CXCR4')], assay = 'RNA', layer = 'data') %>% t

metadata = tibble(
  cellid = names(Idents(objs_by_time$s3299$T3)), 
  assign = unname(Idents(objs_by_time$s3299$T3))
)

df %>% 
  as.data.frame %>% 
  tibble::rownames_to_column('cellid') %>% 
  full_join(., metadata, by = 'cellid') %>% 
  filter(assign %in% c('CXCR4hiCD27lo', 'CXCR4loCD27hi')) %>% 
  # ggplot(aes(x = assign, y = TKT)) + 
  # geom_violin()+
  # geom_jitter() + 
  # geom_point() +
  ggplot(aes(TKT)) + 
  geom_histogram(binwidth = 0.01) + 
  facet_wrap(vars(assign)) + 
  theme_bw()

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

data_dge = list(
  's19' = objs_by_time$s19[c("T4", "T3", "T5", "T1")], 
  's3299' = objs_by_time$s3299[c('T1', 'T2')]
)

dge = lapply(data_dge, function(x){
  lapply(x, function(s) {
    run_dge(s, idents = rev(c('CXCR4hiCD27lo', 'CXCR4loCD27hi')), name_ident = 'annotation_final')
  })
})
saveRDS(dge, 'Hemasphere_data/dge_CXCR4hiCD27lo_vs_CXCR4loCD27hi_s19_s3299.rds')

res_dge_ppp = lapply(dge %>% names, function(x) {
  lapply(dge[[x]] %>% names, function(s) {
    dge[[x]][[s]] %>% 
      tibble::rownames_to_column('gene') %>% 
      mutate(sample = paste(x, s, sep = '_')) %>% 
      filter(gene %in% c(ppp, 'CXCR4', 'CD27'))
  }) %>% bind_rows()
}) %>% bind_rows()

cols = setNames(
  c(kelly.colors(22), 'indianred'), nm = res_dge_ppp$gene %>% unique
)
res_dge_ppp %>% 
  mutate(sign = ifelse(p_val_adj < 0.05, '*', '')) %>% 
  # mutate(sample = factor(sample)) %>% 
  tidyr::separate(sample, into = c('patient', 'time_point'), sep = '_') %>% 
  ggplot(aes(gene, y = avg_log2FC, fill = gene, label = sign)) +
  geom_bar(stat = 'identity', position = 'dodge') + 
  facet_grid(patient ~ time_point, scales = 'free') + 
  theme_light() + 
  scale_fill_manual(values = cols) + 
  ggrepel::geom_text_repel(nudge_x = 0, direction = 'y') +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), 
        legend.position = 'bottom', legend.direction = 'horizontal') + 
  ggtitle('CXCR4hiCD27lo_vs_CXCR4loCD27hi')
ggsave('Hemasphere_data/dge_CXCR4hiCD27lo_vs_CXCR4loCD27hi_s19_s3299.pdf', width = 15, height = 10)

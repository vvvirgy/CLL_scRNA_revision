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

# ppp = c(
#  'ALDOA',	
#  'ALDOB',	
#  'ALDOC',	
#  'DERA',	
#  'FBP1',	
#  # 'FBP2',	
#  'G6PD',	
#  'GPI',	
#  'H6PD',	
#  'PFKL',	
#  'PFKM',	
#  'PFKP',	
#  'PGD',	
#  'PGLS',	
#  'PGM1',	
#  'PGM2',	
#  'PRPS1',	
#  # 'PRPS1L1',	
#  'PRPS2',	
#  'RBKS',	
#  'RPE',	
#  # 'RPEL1',	
#  'RPIA',	
#  'TALDO1',	
#  'TKT',	
#  'TKTL1',	
#  'TKTL2'
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

DotPlot(object = objs_by_time$s19$T3, features = c(ppp, 'CXCR4', 'CD27', 'CD5'), cols = 'RdBu',
        idents = c('CXCR4hiCD27lo', 'CXCR4loCD27hi'), dot.scale = 19) +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1, size = 16), 
        legend.position = 'right', legend.direction = 'vertical', 
        axis.line = element_line(linewidth = 0.3), 
        axis.ticks = element_line(linewidth = 0.3)) +
  xlab('') + 
  ylab('') 
  # guides(color = guide_legend(nrow = 1, byrow = FALSE),
  #        size = guide_legend(nrow = 1, byrow = TRUE)) + 
  # ggtitle('Donor 19-T3')
ggsave('Hemasphere_data/ppp_different_populations_s19_t3.png', width = 13, height = 5)
ggsave('Hemasphere_data/ppp_different_populations_s19_t3.pdf', width = 9, height = 4)

DotPlot(object = objs_by_time$s3299$T1, features = c(ppp, 'CXCR4', 'CD27', 'CD5'), cols = 'RdBu',
        idents = c('CXCR4hiCD27lo', 'CXCR4loCD27hi'), dot.scale = 19) +
  theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1, size = 16), 
        legend.position = 'right', legend.direction = 'vertical', 
        axis.line = element_line(linewidth = 0.3), 
        axis.ticks = element_line(linewidth = 0.3)) +
  xlab('') + 
  ylab('') 
ggsave('Hemasphere_data/ppp_different_populations_s3299_t1.png', width = 10, height = 5)
ggsave('Hemasphere_data/ppp_different_populations_s3299_t1.pdf', width = 9, height = 4)


cells_to_vis = objs_by_time$s19$T3@meta.data %>% 
  dplyr::filter(annotation_final %in% c('CXCR4hiCD27lo', 'CXCR4loCD27hi')) %>% 
  tibble::rownames_to_column('id') %>% 
  group_by(annotation_final) %>% 
  group_split() 
names(cells_to_vis) = lapply(cells_to_vis, function(x) {
  x %>% 
    pull(annotation_final) %>% unique
}) %>% unlist
cells_to_vis = lapply(cells_to_vis, function(x) {x$id})

pdf('Hemasphere_data/s19_t3_umap.pdf')
DimPlot(objs_by_time$s19$T3, cells.highlight = cells_to_vis, cols.highlight = c('CXCR4hiCD27lo' = '#16325B', 'CXCR4loCD27hi' = '#0F828C')) + 
  theme(legend.position = 'bottom')
dev.off()

# visualize umaps 
cells_to_vis = objs_by_time$s3299$T1@meta.data %>% 
  dplyr::filter(annotation_final %in% c('CXCR4hiCD27lo', 'CXCR4loCD27hi')) %>% 
  tibble::rownames_to_column('id') %>% 
  group_by(annotation_final) %>% 
  group_split() 
names(cells_to_vis) = lapply(cells_to_vis, function(x) {
  x %>% 
    pull(annotation_final) %>% unique
}) %>% unlist
cells_to_vis = lapply(cells_to_vis, function(x) {x$id})

pdf('Hemasphere_data/s3299_t1_umap.pdf')
DimPlot(objs_by_time$s3299$T1, cells.highlight = cells_to_vis, cols.highlight = c('CXCR4hiCD27lo' = '#16325B', 'CXCR4loCD27hi' = '#0F828C')) + 
  theme(legend.position = 'bottom')
dev.off()

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

dge_v2 = lapply(data_dge, function(x){
  lapply(x, function(s) {
    run_dge(s, idents = rev(c('CXCR4hiCD27lo', 'CXCR4loCD27hi')), name_ident = 'annotation_final')
  })
})
saveRDS(dge_v2, 'Hemasphere_data/dge_CXCR4hiCD27lo_vs_CXCR4loCD27hi_s19_s3299.rds')

res_dge_ppp = lapply(dge_v2 %>% names, function(x) {
  lapply(dge_v2[[x]] %>% names, function(s) {
    dge_v2[[x]][[s]] %>% 
      tibble::rownames_to_column('gene') %>% 
      mutate(sample = paste(x, s, sep = '_')) %>% 
      filter(gene %in% c(ppp, 'CXCR4', 'CD27', 'CD5'))
  }) %>% bind_rows()
}) %>% bind_rows() 

nadeu_dge = res_dge_ppp %>% 
  dplyr::filter(sample %in% c('s19_T3', 's3299_T1'))

write.table(nadeu_dge, 'Hemasphere_data/ppp_dge_nadeu.csv', sep = ',', quote = F, row.names = F, col.names = T)

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





# test with heatmaps --> singlecellcomplexheatmap?
x$timepoint = x$time_point
x$celltype = x$annotation_final
x$group = paste(x$time_point, x$celltype, sep = '_')

x@meta.data %>% head

gene_groups = list(
  'Pentose phosphate pathway' = ppp, 
  'Cell markers' = c('CD27', 'CXCR4', 'CD5')
)

create_single_cell_complex_heatmap(
  seurat_object = x,
  features = c(ppp, 'CD27', 'CXCR4', 'CD5'),
  group_by = 'group', 
  gene_classification = gene_groups, 
  color_palette = c("navy", "white", "firebrick"), 
  idents = c('CXCR4hiCD27lo', 'CXCR4loCD27hi')
) + 
  theme(legend.position = 'bottom')


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

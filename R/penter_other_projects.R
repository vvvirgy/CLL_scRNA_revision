rm(list=ls())
library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(harmony)
library(patchwork)

setwd('/orfeo/cephfs/scratch/cdslab/vgazziero/TLS/CLL_revision')

good_samples = readRDS('Hemasphere_data/good_samples_penter.rds')

# genes list
genes_list = list(
  'transporters' = c(
    'CEMIP2',
    'SLC1A5',
    'SLC7A5',
    'SCL3A2',
    'SLC38A1',
    'SLC38A2',
    'SLC6A14'), 
  'ILT3' = c(
    'LILRB4'
  ), 
  'copper' = c(
    'SLC31A1', 
    'CD44')
)

# create dotplots 


lapply(good_samples %>% names, function(x) {
  lapply(names(genes_list), function(g) {
    pdf(paste0('other_projects/', g, '_different_populations_', x, '.pdf'), width = 9, height = 5)
    print(DotPlot(object = good_samples[[x]], features = c(genes_list[[g]], 'CXCR4', 'CD27', 'CD5'), cols = 'RdBu',
                  idents = c('CXCR4hiCD27lo', 'CXCR4loCD27hi'), dot.scale = 19) +
            theme(axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1, size = 16), 
                  legend.position = 'right', legend.direction = 'vertical', 
                  axis.line = element_line(linewidth = 0.3), 
                  axis.ticks = element_line(linewidth = 0.3)) +
            xlab('') + 
            ylab('')) 
    dev.off()
  })
})

# perform differential expression analysis
dge = readRDS('Hemasphere_data/dge_penter_all.rds')
genes = unlist(genes_list)

dge_penter = lapply(dge %>% names, function(x) {
  dge[[x]] %>% 
    dplyr::mutate(sample = x) %>% 
    tibble::rownames_to_column('gene') %>% 
    dplyr::filter(gene %in% c(genes, 'CXCR4', 'CD27', 'CD5'))
}) %>% bind_rows()
write.table(dge_penter, 'other_projects/dge_penter.csv', sep = ',', quote = F, row.names = F, col.names = T)

# lapply(good_samples, function(s){
#   
#   cells_to_vis = s@meta.data %>% 
#     dplyr::filter(annotation_final %in% c('CXCR4hiCD27lo', 'CXCR4loCD27hi')) %>% 
#     tibble::rownames_to_column('id') %>% 
#     group_by(annotation_final) %>% 
#     group_split() 
#   
#   names(cells_to_vis) = lapply(cells_to_vis, function(x) {
#     x %>% 
#       pull(annotation_final) %>% unique
#   }) %>% unlist
#   
#   cells_to_vis = lapply(cells_to_vis, function(x) {x$id})
#   
#   sid = s$sample_id %>% unique
#   
#   pdf(paste0('Hemasphere_data/', sid, '_umap.pdf'))
#   print(DimPlot(s, cells.highlight = cells_to_vis, cols.highlight = c('CXCR4hiCD27lo' = '#16325B', 'CXCR4loCD27hi' = '#0F828C')) + 
#           theme(legend.position = 'bottom') + 
#           ggtitle(sid))
#   dev.off()
# })







write.table(ppp_dge_penter, 'Hemasphere_data/ppp_dge_penter.csv', sep = ',', quote = F, row.names = F, col.names = T)

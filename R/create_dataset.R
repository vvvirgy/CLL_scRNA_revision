rm(list=ls())
setwd("/orfeo/cephfs/scratch/cdslab/vgazziero/TLS/CLL_revision")

library("GEOquery")
library(tidyverse)

# data = getGEO(GEO = 'GSE165087', GSEMatrix = FALSE) #, destdir = '/orfeo/LTS/CDSLab/LT_storage/vgazziero/cll_tls_prj')

penter_data = getGEO(GEO = 'GSE165087', GSEMatrix = TRUE, destdir = '/orfeo/LTS/CDSLab/LT_storage/vgazziero/cll_tls_prj') 
metadata = pData(phenoData(penter_data[[1]])) 
metadata$`time:ch1` %>% unique

metadata_v2 <- readxl::read_excel("/orfeo/LTS/CDSLab/LT_storage/vgazziero/cll_tls_prj/Penter_dataset_metadata.xlsx")
metadata_v2 %>% colnames

downloads = metadata_v2 %>% 
  dplyr::filter(Analizzare == 'SI') %>% 
  # pull(supplementary_file_1)
  dplyr::select(title, starts_with('supplementary'), relation) 
write.table(downloads, 'Penter_data/download_samples.csv', sep = ',', quote = F, col.names = T, row.names = F)

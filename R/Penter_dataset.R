rm(list=ls())
library(tidyverse)
library(Seurat)

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
# Visualize QC metrics as a violin plot
VlnPlot(xs$CLL1_3.1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

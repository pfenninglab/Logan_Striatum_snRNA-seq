## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)
library(SingleCellExperiment)
library(DelayedArray)
library(HDF5Array)
library(Matrix.utils)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/HeKleyman2021_macaque_striatum_data_processing'
plan("multicore", workers = 16)
options(future.globals.maxSize = 20 * 1024^3)

#################################################
## 1) load only the scaled and raw RNA counts
h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
if(!file.exists(file.path(h5Dir, 'GSE167920_Results_mergeNuclei_processed_finalassays.h5'))){
obj_all = here(DATADIR, 'rdas', 'GSE167920_Results_full_nuclei_processed_final.h5Seurat') %>% 
  LoadH5Seurat(assay = 'RNA')
obj_msn = here(DATADIR, 'rdas', 'GSE167920_Results_MSNs_processed_final.h5Seurat') %>% 
  LoadH5Seurat(assay = 'RNA')

head(obj_all[[]])
obj_all = obj_all[,!obj_all$cell_type_2 %in% c('MSNs')]
table(obj_all$cell_type_2)

head(obj_msn[[]])
table(obj_msn$MSN_type)

obj = merge(obj_all, y = obj_msn, add.cell.ids = c("Other", "MSN"), 
            project = "Macaque Striatum")

## gather & rename the cell types
obj$celltype3 = ifelse(is.na( obj$MSN_type), obj$cell_type_2,  obj$MSN_type)
obj$celltype3 = ifelse( obj$celltype3 == 'Interneurons', 'Interneuron',   obj$celltype3)
obj$celltype3 = ifelse( obj$celltype3 == 'Mural/Fibroblast', 'Mural',   obj$celltype3)
obj$celltype3 = obj$celltype3 %>% make.names()
table(obj$celltype3)

## convert to SingleCellExperiment
sce = as.SingleCellExperiment(obj)
saveHDF5SummarizedExperiment(sce, h5Dir, prefix="GSE167920_Results_mergeNuclei_processed_final", replace=TRUE)
}

## load hdf5 file w/ the raw counts
sce = loadHDF5SummarizedExperiment(h5Dir, prefix="GSE167920_Results_mergeNuclei_processed_final")

##################################################
# 2) create or load pseudobulk sce object per animal
save_pseudobulk =here(DATADIR, 'rdas', 'GSE167920_Results_mergeNuclei_processed_final.sce.rds')
if(!file.exists(save_pseudobulk)){
  ## aggregate by cluster-sample to create pseudo-bulk count matrix
  colData(sce); table(sce$celltype3)
  groups <- colData(sce)[, c("celltype3", 'split')]
  pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
  dim(pb)

  ## split by cluster, transform & rename columns
  pb_colData = colData(sce) %>% as.data.frame() %>%
    rownames_to_column('match') %>% 
    mutate(match = paste(celltype3, split, sep = '_')) %>% 
    filter(!duplicated(match)) %>% column_to_rownames('match')
  pb_colData = pb_colData[rownames(pb),]
  dim(pb_colData)
  
  ## add number of cells per aggregate
  num_cells = groups %>% as.data.frame() %>% 
    mutate(tmp = paste(celltype3, split,sep= '_')) %>% pull(tmp) %>% table()
  num_cells %>% as.numeric() %>% summary()
  pb_colData$numCells = num_cells[rownames(pb_colData)]
  
  ## add the gene detection rate
  pb_colData$cdr <- scale(rowMeans(pb > 0)) 
  
  ## create SingleCellExperiment from pseudo bulk counts across all cell types and region_name
  (pb <- SingleCellExperiment(assays = t(pb), colData = pb_colData))
  saveRDS(pb, save_pseudobulk)
} else {
  pb = readRDS(save_pseudobulk)
}

table(pb$celltype3)


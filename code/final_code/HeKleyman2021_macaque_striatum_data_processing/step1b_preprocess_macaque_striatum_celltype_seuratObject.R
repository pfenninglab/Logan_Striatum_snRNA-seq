## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/HeKleyman2021_macaque_striatum_data_processing'
plan("multicore", workers = 8)
options(future.globals.maxSize = 20000 * 1024^2)

#######################################################################
# 1) read in labeled monkey all nuclei dataset, He, Kleyman et al. 2021
## read in the all nuclei object
monkey_all = here(DATADIR, 'rdas', 'GSE167920_Results_full_nuclei_processed_final.rds') %>% 
  readRDS() 
monkey_all$split = with(monkey_all[[]], 
                        paste(monkey, region_name, gsub('[ACGT-]','', colnames(monkey_all)),sep = '-'))

## split by monkey, brain region, and compute SCTransform
monkey_all_list <- SplitObject(monkey_all, split.by = "split")
monkey_all_list <- lapply(X = monkey_all_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = monkey_all_list, nfeatures = 3000)
monkey_all_list <- PrepSCTIntegration(object.list = monkey_all_list, anchor.features = features)
monkey_all_list <- lapply(X = monkey_all_list, FUN = RunPCA, features = features)

## co-embed data into 1 SCTransformed space
monkey_all.anchors <- FindIntegrationAnchors(
  object.list = monkey_all_list, normalization.method = "SCT",
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
monkey_all.combined.sct <- IntegrateData(anchorset = monkey_all.anchors, normalization.method = "SCT", dims = 1:30)
monkey_all.combined.sct <- RunPCA(monkey_all.combined.sct, verbose = FALSE)
monkey_all.combined.sct <- RunUMAP(monkey_all.combined.sct, reduction = "pca", dims = 1:30)

## save SCTransformed all counts matrix
save_all_h5 = here(DATADIR, 'rdas', 'GSE167920_Results_full_nuclei_processed_final.h5Seurat')
SaveH5Seurat(monkey_all.combined.sct, filename = save_all_h5, overwrite = TRUE)

#######################################################
# 2) read in labeled monkey MSN subtype dataset, He, Kleyman et al. 2021
## read in the MSN sub-type  object
monkey_msn = here(DATADIR, 'rdas', 'GSE167920_Results_MSNs_processed_final.rds') %>% 
  readRDS()
monkey_msn$split = with(monkey_msn[[]], paste(monkey, region_name))

## split by monkey, brain region, and comput SCTransform
monkey_msn_list <- SplitObject(monkey_msn, split.by = "split")
monkey_msn_list <- lapply(X = monkey_msn_list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = monkey_msn_list, nfeatures = 3000)
monkey_msn_list <- PrepSCTIntegration(object.list = monkey_msn_list, anchor.features = features)
monkey_msn_list <- lapply(X = monkey_msn_list, FUN = RunPCA, features = features)

## co-embed data into 1 SCTransformed space
monkey_msn.anchors <- FindIntegrationAnchors(
  object.list = monkey_msn_list, normalization.method = "SCT",
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
monkey_msn.combined.sct <- IntegrateData(anchorset = monkey_msn.anchors, normalization.method = "SCT", dims = 1:30)
monkey_msn.combined.sct <- RunPCA(monkey_msn.combined.sct, verbose = FALSE)
monkey_msn.combined.sct <- RunUMAP(monkey_msn.combined.sct, reduction = "pca", dims = 1:30)

## save SCTransformed MSN counts matrix
save_msn_h5 = here(DATADIR, 'rdas', 'GSE167920_Results_MSNs_processed_final.h5Seurat')
SaveH5Seurat(monkey_msn.combined.sct, filename = save_msn_h5, overwrite = TRUE)

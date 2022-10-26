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

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 28 cores
plan("multicore", workers = 28)
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

##################################################
# 1) load in cell type labels for label transfer
## read in Logan BU snRNA dataset to label transfer
save_merged_fn = here('data/tidy_data/Seurat_projects', 
                        "OUD_Striatum_filtered_SCT_SeuratObj_N22.h5Seurat")
obj_merged = save_merged_fn %>% LoadH5Seurat() 

## read in labeled monkey MSN subtype dataset, He, Kleyman et al. 2021
## subject to same region in monkey, recompute PCA space
monkey_all = here(DATADIR, 'rdas', 
                  'GSE167920_Results_full_nuclei_processed_final.h5Seurat') %>% 
  LoadH5Seurat(assays = "integrated") %>% 
  subset(subset = region_name %in% c( "caudate", 'putamen')) %>% 
  RunPCA(verbose = FALSE)

head(monkey_all[[]])

## dorsal striatum subtypes
subtypes = c('D1-Matrix', 'D2-Matrix', 'D1-Striosome', 'D2-Striosome', 'D1/D2-Hybrid')
monkey_msn = here(DATADIR, 'rdas', 'GSE167920_Results_MSNs_processed_final.h5Seurat') %>% 
  LoadH5Seurat(assays = "integrated") %>% 
  subset(subset = region_name %in% c( "caudate", 'putamen')) %>% 
  subset(subset = MSN_type %in% subtypes) %>%
  RunPCA(verbose = FALSE)

head(monkey_msn[[]])


###############################################################
# 2) transfer cell type labels from Macaque all nuclei dataset

## compute anchor cells b/t BU dataset and Macaque caudate dataset
anchors.mac_all <- FindTransferAnchors(
  reference = monkey_all, query = obj_merged, reduction = 'rpca',
  query.assay = 'integrated', reference.assay = 'integrated',
  normalization.method = 'SCT', dims = 1:30, reference.reduction = "pca")

## predict BU cell type w/ macaque 'cell_type_2' column
predictions.mac_all <- TransferData(
  anchorset = anchors.mac_all, refdata = monkey_all$cell_type_2, dims = 1:30, store.weights = F)
score_threshold = .5
predictions.mac_all$predicted.id = 
  with(predictions.mac_all, ifelse(prediction.score.max < score_threshold, 'UNK_ALL', predicted.id))

## append the predicted cell type to 
obj_merged <- AddMetaData(obj_merged, metadata = predictions.mac_all$predicted.id,
                            col.name = 'celltype1')

## see which Seurat clusters corresponds to macaque cell types
(tbl1 = with(obj_merged[[]], table(celltype1, seurat_clusters)))
cluster_to_macAll = colnames(tbl1)
names(cluster_to_macAll) = rownames(tbl1)[apply(tbl1, 2, which.max)]
cluster_to_macAll


######################################################################
# 3) transfer cell type labels from Macaque MSN sub-type nuclei dataset
## split the BU data to just predicted MSN subtypes
## do this with pipes in one go
obj_msn = obj_merged %>% 
  subset(subset = celltype1 == 'MSNs' | 
           seurat_clusters %in% cluster_to_macAll[names(cluster_to_macAll) =='MSNs']) %>% 
  RunPCA(verbose = FALSE) %>% 
  FindNeighbors(dims = 1:30, verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, algorithm = 2,verbose = FALSE)

## reciprocal PCA project b/t monkey/human MSNs
anchors.mac_msn <- FindTransferAnchors(
  reference = monkey_msn, query = obj_msn, reduction = 'rpca', 
  query.assay = 'integrated', reference.assay = 'integrated',
  normalization.method = 'SCT',  dims = 1:30, reference.reduction = "pca")

## predict BU MSNs at the subtype level
predictions.mac_msn <- TransferData(anchorset = anchors.mac_msn, 
                                    refdata = monkey_msn$MSN_type, dims = 1:30)

## set less confident predictions to unknown celltype
score_threshold = 0.5
predictions.mac_msn$predicted.id =
  with(predictions.mac_msn, ifelse(prediction.score.max < score_threshold, 'UNK_MSN', predicted.id))

## transfer MSN subtypes labels back to MSN and full dataset
obj_msn <- AddMetaData(obj_msn, metadata = predictions.mac_msn$predicted.id,
                            col.name = 'celltype2')
obj_merged$celltype2 = obj_msn[['celltype2']][colnames(obj_merged),]
obj_merged$celltype2 = with(obj_merged[[]], ifelse(is.na(celltype2 ),celltype1, celltype2 ))

## see which cluster clusters correponds to MSNs
(tbl2 = with(obj_msn[[]], table(celltype2, seurat_clusters)))
cluster_to_macMSN = colnames(tbl2)
names(cluster_to_macMSN) = rownames(tbl2)[apply(tbl2, 2, which.max)]
table(obj_msn$celltype2)
cluster_to_macMSN


############################################
# 4) save the subset of cells that are MSNs
save_subset_msn = here('data/tidy_data/Seurat_projects',
                "OUD_Striatum_subsetMSN_SCT_SeuratObj_N22.h5Seurat")
SaveH5Seurat(obj_msn, filename = save_subset_msn, overwrite = TRUE)
SaveH5Seurat(obj_merged, filename = save_merged_fn, overwrite = TRUE)


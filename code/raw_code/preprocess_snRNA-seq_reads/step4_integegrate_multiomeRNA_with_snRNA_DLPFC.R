library(scds)
library(SingleCellExperiment)
library(DropletUtils)
library(Seurat)
library(Matrix)
library(tidyverse)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

############################################################
## load SeuratObjects of labeled all nuclei, from snRNA-seq
DATADIR1='data/tidy_data/rdas/JH_labeled_PFC_LabeledNuclei_202120827'
obj_rna = readRDS(file.path(DATADIR1, 'nuclei_all_final.rds'))
obj_rna_exc = readRDS(file.path(DATADIR1, 'excita_final.rds'))
obj_rna_inh = readRDS(file.path(DATADIR1, 'interneurons_final.rds'))


#############################################
## transfer neuronal subtypes to full dataset
obj_rna@meta.data$cell_type = as.character(obj_rna@meta.data$cell_class)
indEXC = match(rownames(obj_rna_exc@meta.data), rownames(obj_rna@meta.data))
indINH = match(rownames(obj_rna_inh@meta.data), rownames(obj_rna@meta.data))

obj_rna@meta.data$cell_type[indEXC] = as.character(obj_rna_exc@meta.data$cell_type)
obj_rna@meta.data$cell_type[indINH] = as.character(obj_rna_inh@meta.data$cell_type) 
obj_rna@meta.data$cell_type = gsub(' [1-9]$', '', obj_rna@meta.data$cell_type)
table(obj_rna@meta.data$cell_type)

# some cell types not labeled in sub analyses, drop
obj_rna@meta.data$cell_type[obj_rna@meta.data$cell_type %in% c('Excita_neurons', 'Interneurons')] = 'Drop'
obj_rna = subset(x = obj_rna, subset = cell_type != 'Drop')


#############################################
## load in raw counts of unlabeled multiomeRNA
sce_fn = file.path('data/tidy_data/rdas', 'raw_multiomeRNA_DLPFC_sce.rds')
sce = readRDS(sce_fn)
rownames(sce) = gsub('gene-', '', rownames(rowData(sce)))
rownames(sce) = gsub('-[1-9]$', '', rownames(rowData(sce)))
sce = sce[!duplicated(rownames(sce))]

obj_multi = as.Seurat(sce, data = 'counts')
obj_multi <- NormalizeData(obj_multi, verbose = FALSE)
obj_multi <- FindVariableFeatures(obj_multi, selection.method = "vst", 
                                  nfeatures = 2000, verbose = FALSE)


#############################################
## rPCA integration of snRNA and multiomeRNA
obj_list = list('snRNA'= obj_rna, 'multiome' = obj_multi)

features <- SelectIntegrationFeatures(object.list = obj_list)
obj_list <- lapply(X = obj_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# find integration anchors
obj_anchors <- FindIntegrationAnchors(
  object.list = obj_list, anchor.features = features, 
  reduction = "rpca") #, k.anchor = 20

# this command creates an 'integrated' data assay
obj_combined <- IntegrateData(anchorset = obj_anchors)
obj_combined <- ScaleData(obj_combined, verbose = FALSE)
obj_combined <- RunPCA(obj_combined, npcs = 30, verbose = FALSE)
obj_combined <- RunUMAP(obj_combined, reduction = "pca", dims = 1:30)
obj_combined <- FindNeighbors(obj_combined, reduction = "pca", dims = 1:30)
obj_combined <- FindClusters(obj_combined, resolution = 0.5)

integrated_fn = file.path('data/tidy_data/rdas', 
                          'snRNA_multiomeRNA_integrated_macaque_DLPFC_N60k.rds')
saveRDS(obj_combined, integrated_fn)


############################################################################
## transfer labels from the monkey O & S snRNA dataset to multiomeRNA dataset
trans_anchors <- FindTransferAnchors(reference.reduction = "pca", 
  reference = obj_rna, query = obj_multi, dims = 1:30)
predictions <- TransferData(anchorset = trans_anchors, 
                            refdata = obj_rna$cell_type, dims = 1:30)

obj_multi <- AddMetaData(obj_multi, metadata = predictions)
transferred_fn = file.path('data/tidy_data/rdas', 
                          'multiomeRNA_transferredLabels_macaque_DLPFC_N30k.rds')
saveRDS(obj_multi, transferred_fn)

metadata_fn = file.path('data/tidy_data/rdas', 
                           'multiomeRNA_transferredLabels_macaque_DLPFC_N30k.meta.rds')
meta = obj_multi@meta.data
rownames(meta) = with(meta, paste0(Sample, '#', Barcode))
saveRDS(meta, metadata_fn)




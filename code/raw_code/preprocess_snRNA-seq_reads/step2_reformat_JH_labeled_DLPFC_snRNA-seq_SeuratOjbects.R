# conda activate r4
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(Seurat)
library(SingleCellExperiment)
library(tidyverse)
library(Matrix)

DATADIR='data/tidy_data/rdas/JH_labeled_PFC_LabeledNuclei_202120827'

## load SeuratObjects of labeled excitatory and inhibitory neurons
obj_exc = readRDS(file.path(DATADIR, 'excita_final.rds'))
obj_inh = readRDS(file.path(DATADIR, 'interneurons_final.rds'))
obj_neuron <- merge(obj_exc, y = obj_inh, add.cell.ids = c("EXC", "INH"), 
                    project = "DLPFC_snRNA")

table(obj_neuron@meta.data$cell_type)

sce_neuron = as.SingleCellExperiment(obj_neuron, assay = 'RNA')
saveRDS(sce_neuron, file.path(DATADIR, 'neuron_final.sce.rds'))

## load SeuratObjects of labeled all nuclei, glia
obj = readRDS(file.path(DATADIR, 'nuclei_all_final.rds'))
table(obj@meta.data$cell_class)

gliaTypes = c('Astrocytes', 'Oligos','Microglia', 'Oligo_Pre', 'Endothelial')
obj_glia = subset(obj, subset = cell_class %in% gliaTypes )

sce_glia = as.SingleCellExperiment(obj_glia, assay = 'RNA')
saveRDS(sce_glia, file.path(DATADIR, 'glia_final.sce.rds'))

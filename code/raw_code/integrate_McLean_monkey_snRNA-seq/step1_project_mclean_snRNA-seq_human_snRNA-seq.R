## conda activate r4
## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## packages for data table processing 
library(here)
library(tidyverse)
library(writexl)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 28 cores
plan("multicore", workers = 12)
options(future.globals.maxSize = 80 * 1024^3)

neuron_types = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid', 
                 'D1-ICj', 'D1-NUDAP', 'D1-Shell/OT', 'D2-Shell/OT',
                 'Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH')
glia_types = c( 'Astrocytes',  'Microglia', 'Oligos', 'Oligos_Pre')

##############################################
# 1) import the seurat objects to be projected 
# read in Logan snRNA dataset to plot
obj_human = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat() 
obj_human = obj_human %>% RunUMAP(dims = 1:30, return.model = TRUE)

obj_hg_neuron = obj_human %>% 
  subset(celltype3 %in% neuron_types) %>% 
  RunUMAP(dims = 1:30, return.model = TRUE)

obj_hg_glia = obj_human %>% 
  subset(celltype3 %in% glia_types) %>% 
  RunUMAP(dims = 1:30, return.model = TRUE)

# read in McLean snRNA dataset to plot
obj_monkey = file.path('/projects/pfenninggroup/singleCell/McLean_chronic_opioid_monkey_snRNA-seq',
                       'data/tidy_data/Seurat_projects/OUD_Striatum_refined_all_SeuratObj_N16.h5Seurat') %>% 
  LoadH5Seurat() 

obj_rm_neuron = obj_monkey %>% subset(celltype3 %in% neuron_types)
obj_rm_glia = obj_monkey %>% subset(celltype3 %in% glia_types)

table(obj_monkey$celltype3)


########################################################################
# 2) project the monkey neuronal and glial snRNA-seq into the human UMAP
anchors.neuron <- FindTransferAnchors(reference = obj_hg_neuron, query = obj_rm_neuron,
                               dims = 1:30, reference.reduction = "pca")

query.neuron <- MapQuery(anchorset = anchors.neuron, reference = obj_hg_neuron, 
                  query = obj_rm_neuron, refdata = list(celltype = "celltype3"), 
                  reference.reduction = "pca", reduction.model = "umap")
table(query.neuron$celltype3, query.neuron$predicted.celltype)



anchors.glia <- FindTransferAnchors(reference = obj_hg_glia, query = obj_rm_glia,
                                      dims = 1:30, reference.reduction = "pca")

query.glia <- MapQuery(anchorset = anchors.glia, reference = obj_hg_glia, 
                         query = obj_rm_glia, refdata = list(celltype = "celltype3"), 
                         reference.reduction = "pca", reduction.model = "umap")
table(query.glia$celltype3, query.glia$predicted.celltype)



p1 <- DimPlot(obj_human, reduction = "umap", group.by = "celltype3", label = TRUE, 
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Human annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query monkey annotations")
p1 + p2

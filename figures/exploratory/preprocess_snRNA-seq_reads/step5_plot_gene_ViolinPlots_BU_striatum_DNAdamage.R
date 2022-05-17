## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## gene expression heatmap plots
library(Nebulosa)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

PLOTDIR='figures/exploratory/preprocess_snRNA-seq_reads'

###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid', 'UNK_MSN')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a', '#b15928')
names(subtypes_col) = subtypes

othertypes = c('Interneurons', 'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural/Fibroblast', 'Oligos', 'Oligos_Pre', 'UNK_ALL')
othertypes_col = carto_pal(length(othertypes), "Vivid")
names(othertypes_col) = othertypes


##################################################
# 1) read in per-cell DNA-damage scores
df = readRDS(here('data/tidy_data/Seurat_projects/AUCell/22_02_04_cells_rankings.rds'))



## plot the opioid receptors and ligands
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_all.DNAdamage.pdf'), width = 7.25, height = 12)
VlnPlot(obj_merged, features =c(markDNAdamage), ncol = 3, slot = "data", pt.size = 0,
        fill.by = 'celltype1', group.by = 'celltype1', 
        cols = c('MSNs' = 'gray', othertypes_col), log = T) & 
  theme(legend.position = 'none', axis.title=element_blank()) & 
  coord_flip() 
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_all.DNAdamageByDxSUD.pdf'), width = 7.25, height = 12)
VlnPlot(obj_merged, features =c(markDNAdamage), slot = "data", pt.size = 0, 
        ncol = 3, group.by = 'celltype1', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), 
        split.plot = TRUE, log = T) & 
  theme(axis.title=element_blank()) & coord_flip() 
dev.off()


#####################################
# 2) plot genes at MSN subtype level
## read in Logan BU snRNA dataset to label transfer
obj_msn = here('data/tidy_data/Seurat_projects', 
               "BU_Run1_Striatum_subsetMSN_SCT_SeuratObj_N4.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')
obj_msn$celltype2 = factor(obj_msn$celltype2 , c(subtypes, othertypes))
obj_msn = obj_msn[, !obj_msn$celltype2 %in% c('UNK_MSN', 'UNK_ALL')]
Idents(obj_msn) = 'celltype2'

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_msn.DNAdamage.pdf'), width = 9, height = 12)
VlnPlot(obj_msn, features =c(markDNAdamage), slot = "data", cols = subtypes_col, pt.size = 0) & 
  theme(legend.position = 'none', axis.title=element_blank()) & coord_flip() 
dev.off()


pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_msn.DNAdamageByDxSUD.pdf'), width = 9, height = 12)
VlnPlot(obj_msn, features =c(markDNAdamage), slot = "data", ncol = 3, pt.size = 0,
        split.by = 'DSM.IV.SUD', cols = c('gray', 'red'), split.plot = TRUE) & 
  theme(legend.position = 'none', axis.title=element_blank()) & coord_flip() 
dev.off()



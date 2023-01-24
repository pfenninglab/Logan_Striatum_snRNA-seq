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
# 1) load in full dataset cell type labels for plot

## read in Logan BU snRNA dataset to label transfer
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_filtered_SCT_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')
obj_merged[['group']] = with(obj_merged[[]], ifelse(celltype2 %in% subtypes, 'MSN', 'Other'))
obj_merged = obj_merged[, !obj_merged$celltype1 %in% c('UNK_MSN', 'UNK_ALL')]
obj_merged$group = relevel(factor(obj_merged$group), ref = 'Other')
obj_merged$celltype1 = factor(obj_merged$celltype1 , c('MSNs', othertypes))
obj_merged$celltype2 = factor(obj_merged$celltype2 , c(subtypes, othertypes))
Idents(obj_merged) = 'celltype1'

##################################################
# 2) plot some interesting genes, VGLUT1 and TH
## plot the TH_SLC17A6_ genes
markTH_SLC17A6_ = c('TH', 'SLC17A6')
obj_merged$celltype1 = factor(obj_merged$celltype1, levels = c('MSNs', othertypes)) 
Idents(obj_merged) = 'celltype1'

pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_Viol_all.TH_SLC17A6.pdf'), height = 7.25, width = 12)
VlnPlot(obj_merged, features =c(markTH_SLC17A6_), ncol = 1, slot = "data", pt.size = 0,
        fill.by = 'celltype1', group.by = 'celltype1', 
        cols = c('MSNs' = 'gray', othertypes_col), log = F) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
dev.off()

pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_Viol_all.TH_SLC17A6_ByDxSUD.pdf'), height =7.25,  width = 12)
VlnPlot(obj_merged, features =c(markTH_SLC17A6_), slot = "data", pt.size = 0, 
        ncol = 1, group.by = 'celltype1', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), 
        split.plot = TRUE, log = T) & 
  theme(axis.title=element_blank())
dev.off()



##################################################
# 3) plot some interesting genes, APOE and MTRNR2L8
mark_APOE_MTRNR2L8 = c('APOE', 'MTRNR2L8')

pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_Viol_all.APOE_MTRNR2L8.pdf'), height = 7.25, width = 12)
VlnPlot(obj_merged, features =c(mark_APOE_MTRNR2L8), ncol = 1, slot = "data", pt.size = 0,
        fill.by = 'celltype1', group.by = 'celltype1', 
        cols = c('MSNs' = 'gray', othertypes_col), log = T) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
dev.off()

pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_Viol_all.APOE_MTRNR2L8_ByDxSUD.pdf'), height =7.25,  width = 12)
VlnPlot(obj_merged, features =c(mark_APOE_MTRNR2L8), slot = "data", pt.size = 0, 
        ncol = 1, group.by = 'celltype1', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), 
        split.plot = TRUE, log = T) & 
  theme(axis.title=element_blank())
dev.off()


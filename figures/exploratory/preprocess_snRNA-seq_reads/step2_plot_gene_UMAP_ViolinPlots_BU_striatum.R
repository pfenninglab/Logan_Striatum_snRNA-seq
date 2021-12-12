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
                  "BU_Run1_Striatum_filtered_SCT_SeuratObj_N4.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')
obj_merged[['group']] = with(obj_merged[[]], ifelse(celltype2 %in% subtypes, 'MSN', 'Other'))
obj_merged$group = relevel(factor(obj_merged$group), ref = 'Other')
obj_merged$celltype1 = factor(obj_merged$celltype1 , c('MSNs', othertypes))
obj_merged$celltype2 = factor(obj_merged$celltype2 , c(subtypes, othertypes))
Idents(obj_merged) = 'celltype1'

## Neuron vs. Glia markers
markerGenes1 <- c('RBFOX3', 'LHX6', 'AQP4', 'CX3CR1', 'PDGFRA', 'MOG' )
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_all.allMarkers.pdf'), width = 7.25, height = 4)
p1 = plot_density(obj_merged, features = markerGenes1,  reduction = "umap") +
  plot_layout(nrow = 2, guides = "auto") &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 10))
p1
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_all.allMarkers.pdf'), width = 7.25, height = 5)
VlnPlot(obj_merged, features =c(markerGenes1), slot = "data",  pt.size = 0,
        cols = c('MSNs' = 'gray', othertypes_col)) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
dev.off()


## MSN subtypes
markMSN1 = c('Drd1','Tac1','Reln') %>% toupper()# D1 markers
markMSN2 = c('Drd2','Adora2a','Penk')%>% toupper() # D2 markers
markMSN3 = c('Foxp2', 'Rxfp1', 'Casz1')%>% toupper() # D1/2H markers

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_all.MSNclusterMarkers.pdf'), width = 7.25, height = 3)
p2 = plot_density(obj_merged, slot = 'data', features = markMSN1,  reduction = "umap", joint= T) &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 8))
p3 = plot_density(obj_merged, slot = 'data', features = markMSN2,  reduction = "umap", joint= T) &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 8))
p4 = plot_density(obj_merged, slot = 'data', features = markMSN3,  reduction = "umap", joint= T) &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 8))
p2 + plot_layout(nrow = 1) & theme(legend.position = 'bottom')
p3 + plot_layout(nrow = 1) & theme(legend.position = 'bottom')
p4 + plot_layout(nrow = 1) & theme(legend.position = 'bottom')
dev.off()


## MSN compartments
markMSN4 = c('STXBP6', 'SEMA3E', 'EPHA4', 'GDA') # Matrix markers
markMSN5 =c( 'PDYN', 'OPRM1', 'KHDRBS3', 'KCNIP1') # Striosome markers

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_all.MSNcompMarkers.pdf'), width = 7.25, height = 3)
p5 = plot_density(obj_merged, slot = 'data', features = markMSN4,  reduction = "umap") 
p6 = plot_density(obj_merged, slot = 'data', features = markMSN5,  reduction = "umap", pal = 'inferno')
p5 + plot_layout(nrow = 1) & theme_classic(base_size = 7) & 
  theme(plot.title = element_text(size = 8)) & theme(legend.position = 'bottom') 
p6 + plot_layout(nrow = 1) & theme_classic(base_size = 7) & 
  theme(plot.title = element_text(size = 8)) & theme(legend.position = 'bottom') 
dev.off()


## plot the opioid receptors and ligands
markOPR1 = c('OPRD1', 'OPRM1', 'OPRK1') # opioid receptor genes
markOPR2 = c( 'PENK', 'PDYN') # endorphin peptide genes
obj_merged$celltype1 = factor(obj_merged$celltype1, levels = c('MSNs', othertypes)) 
Idents(obj_merged) = 'celltype1'

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_all.opioids.pdf'), width = 7.25, height = 5)
VlnPlot(obj_merged, features =c(markOPR1, markOPR2), ncol = 3, slot = "data", pt.size = 0,
        fill.by = 'celltype1', group.by = 'celltype1', 
        cols = c('MSNs' = 'gray', othertypes_col), log = T) & 
  theme(legend.position = 'none', axis.title=element_blank()) & 
  coord_flip() 
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_all.opioidsByDxSUD.pdf'), width = 7.25, height = 5)
VlnPlot(obj_merged, features =c(markOPR1, markOPR2), slot = "data", pt.size = 0, 
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
Idents(obj_msn) = 'celltype2'


pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_msn.MSNclusterMarkers.pdf'), width = 7.25, height = 2.5)
VlnPlot(obj_msn, features =c(markMSN1), slot = "data", ncol = 3, pt.size = 0, 
        group.by = 'celltype2', cols = c(subtypes_col, othertypes_col)) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
VlnPlot(obj_msn, features =c(markMSN2), slot = "data", ncol = 3, pt.size = 0, 
        group.by = 'celltype2', cols = c(subtypes_col, othertypes_col)) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
VlnPlot(obj_msn, features =c(markMSN3), slot = "data", ncol = 3, pt.size = 0,
        group.by = 'celltype2', cols = c(subtypes_col, othertypes_col)) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
dev.off()


pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_msn.MSNcompMarkers.pdf'), width = 7.25, height = 5)
VlnPlot(obj_msn, features =c(markMSN4, markMSN5), slot = "data", ncol = 4, pt.size = 0, 
        group.by = 'celltype2', cols = c(subtypes_col, othertypes_col)) & 
  theme(legend.position = 'none', axis.title=element_blank())
dev.off()


pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_msn.opioids.pdf'), width = 7.25, height = 5)
VlnPlot(obj_msn, features =c(markOPR1, markOPR2), slot = "data", cols = subtypes_col, pt.size = 0) & 
  theme(legend.position = 'none', axis.title=element_blank()) & coord_flip() 
dev.off()


pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_Viol_msn.opioidsByDxSUD.pdf'), width = 7.25, height = 5)
VlnPlot(obj_msn, features =c(markOPR1,markOPR2), slot = "data", ncol = 3, pt.size = 0,
        split.by = 'DSM.IV.SUD', cols = c('gray', 'red'), split.plot = TRUE) & 
  theme(legend.position = 'none', axis.title=element_blank()) & coord_flip() 
dev.off()



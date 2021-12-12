## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)
## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

PLOTDIR='figures/exploratory/preprocess_snRNA-seq_reads'

###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix','D1-Striosome', 'D2-Striosome', 'D1/D2-Hybrid', 'UNK_MSN')
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

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_all.ident.pdf'), width = 7.25, height = 4)
Idents(obj_merged) = obj_merged$orig.ident
DimPlot(object = obj_merged, reduction = "umap", split.by = 'orig.ident', 
             cols = brewer.pal(n = 4, name = 'Set2'), label.size = 3 ) +
  guides(color = guide_legend(nrow = 1, override.aes= list(size = 2))) + 
  theme(legend.position = 'bottom')
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_all.clusters.pdf'), width = 7.25, height = 4)
DimPlot(object = obj_merged, reduction = "umap", label.size = 3, 
        group.by = 'seurat_clusters',  split.by = 'orig.ident') +
  guides(color = guide_legend(nrow = 2, override.aes= list(size = 2)))+ 
  theme(legend.position = 'bottom', legend.spacing.x = unit(1.0, 'cm'))
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_all.macaqueLabels.pdf'), width = 7.25, height = 4)
DimPlot(object = obj_merged, reduction = "umap",  split.by = 'orig.ident', 
        group.by = 'celltype2', cols = c(subtypes_col, othertypes_col), label.size = 3) +
  guides(color = guide_legend(nrow = 3, override.aes= list(size = 2))) + 
  theme(legend.position = 'bottom')
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_all.macaqueLabels2.pdf'), width = 7.25, height = 4)
DimPlot(object = obj_merged, reduction = "umap", group.by = 'celltype2', 
        shape.by = 'celltype2', pt.size = .75, label.size = 3) +
  scale_colour_manual(name = "Cell type", labels = names(c(subtypes_col, othertypes_col)),
                      values = c(subtypes_col, othertypes_col)) +   
  scale_shape_manual(name = "Cell type", labels = names(c(subtypes_col, othertypes_col)),
                     values = c(rep(17, 6), rep(19, 8)))
dev.off()


#####################################
# 2) load in MSN sub type labels plot
## read in Logan BU snRNA dataset to label transfer
obj_msn = here('data/tidy_data/Seurat_projects', 
                    "BU_Run1_Striatum_subsetMSN_SCT_SeuratObj_N4.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

## plot MSN subtypes by subjects
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_msn.ident.pdf'), width = 7.25, height = 4)
DimPlot(object = obj_msn, reduction = "umap", group.by = 'orig.ident', label.size = 3, 
        split.by = 'orig.ident', cols = brewer.pal(n = 4, name = 'Set2') ) +
  guides(color = guide_legend(nrow = 1, override.aes= list(size = 2))) + 
  theme(legend.position = 'bottom')
dev.off()

## plot MSN subtypes by Seurat clusters
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_msn.clusters.pdf'), width = 7.25, height = 4)
DimPlot(object = obj_msn, reduction = "umap",  group.by = 'seurat_clusters',
        split.by = 'orig.ident', label.size = 3) +
  guides(color = guide_legend(nrow = 1, override.aes= list(size = 2))) + 
  theme(legend.position = 'bottom')
dev.off()

## plot MSN subtypes by Macaque label transfer
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_msn.macaqueLabels.pdf'), width = 7.25, height = 4)
DimPlot(object = obj_msn, reduction = "umap", split.by = 'orig.ident', 
        group.by = 'celltype2',  cols = subtypes_col, label.size = 3) +
  guides(color = guide_legend(nrow = 2, override.aes= list(size = 2))) + 
  theme(legend.position = 'bottom')
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_UMAP_msn.macaqueLabels2.pdf'), width = 7.25, height = 4)
DimPlot(object = obj_msn, reduction = "umap", group.by = 'celltype2', 
        shape.by = 'celltype2', pt.size = 1, label.size = 3) +
  scale_colour_manual(name = "Cell type", labels = names(c(subtypes_col)),
                      values = c(subtypes_col)) +   
  guides(color = guide_legend(ncol = 1, override.aes= list(size = 2))) + 
  scale_shape_manual(name = "Cell type", labels = names(c(subtypes_col)),
                     values = c(rep(17, 6)))
dev.off()


## plot number of estimated miQC compromised cells
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_CellTypeBarChart.perSample.pdf'), 
    width = 7.25, height = 4, onefile = F)
p1 = obj_merged[[]] %>% ggplot(aes(x = orig.ident, fill = celltype1)) + 
  geom_bar(stat = 'count', position = 'fill')  + ylab('Proportion of cells') + 
  scale_fill_manual(values = c('MSNs' = 'black', subtypes_col, othertypes_col), name = 'Cell type') +
  theme_classic() + xlab('Sample')
p2 =  obj_msn[[]] %>% ggplot(aes(x = orig.ident, fill = celltype2)) + 
  geom_bar(stat = 'count', position = 'fill') + 
  scale_fill_manual(values = c('MSNs' = 'black', subtypes_col, othertypes_col), name = 'Cell type') + 
  theme_classic() + ylab('Proportion of cells') + 
  xlab('Sample')

ggarrange(p1, p2, labels = c("A", "B"), 
          common.legend = TRUE, nrow = 1, legend="bottom")
dev.off()



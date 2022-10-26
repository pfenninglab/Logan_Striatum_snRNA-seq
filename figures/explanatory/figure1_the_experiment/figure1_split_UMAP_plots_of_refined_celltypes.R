ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(ggsci)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

library(here)

DATADIR='data/raw_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure1_the_experiment'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = c(carto_pal(4, "Safe"), 
                   carto_pal(length(othertypes) -4 , "Vivid"))
names(othertypes_col) = othertypes
typecolors = c(subtypes_col, othertypes_col)

## plotting aesthetics
in2mm<-25.4
my_theme = theme_classic(base_size = 6)

###################################
# 1) read in Logan snRNA dataset to plot
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA') 
names(obj_merged[[]] )
table(obj_merged$celltype3 )

fig1_split_umap_allCells_fn = 
  here(PLOTDIR, 'plots', 'fig1_split_umap_allCells.pdf')

pdf(fig1_split_umap_allCells_fn, width = 140/in2mm, height =  100/in2mm)
DimPlot(object = obj_merged, reduction = "umap", split.by = 'DSM.IV.OUD', 
        group.by = 'celltype3', label = T,
        label.size = 1.8, pt.size = 2, ncol = 6, cols = typecolors, raster = T) +
  my_theme +
  guides(color = guide_legend(override.aes = list(size = 2), byrow = F)) +
  theme(legend.position = 'bottom', plot.title= element_blank(),
        legend.spacing.x = unit(6, 'mm'),
        legend.key.size = unit(3, "mm"))
dev.off()



###################################
# 2) plot the avg. cell-wise QC metrics
obj_msn = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_msn_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA') 
names(obj_msn[[]] )
table(obj_msn$celltype3 )

fig1_split_umap_MSNs_fn = 
  here(PLOTDIR, 'plots', 'fig1_split_umap_MSNs.pdf')

pdf(fig1_split_umap_MSNs_fn, width = 140/in2mm, height =  50/in2mm)
DimPlot(object = obj_msn, reduction = "umap", split.by = 'DSM.IV.OUD', 
        group.by = 'celltype3', label = T,
        label.size = 1.5, pt.size = 2, ncol = 6, cols = subtypes_col, raster = T) +
  my_theme +
  theme(legend.position = 'none', plot.title= element_blank())
dev.off()



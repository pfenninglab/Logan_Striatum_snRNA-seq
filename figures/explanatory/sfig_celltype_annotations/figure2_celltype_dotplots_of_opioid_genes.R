ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(ggsci)
library(grid)
library(gridExtra) 

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## num parallel threads for Seurat
plan(sequential) # no parallel
options(future.globals.maxSize = 80e9)
library(here)

DATADIR='data/raw_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure2_the_celltype_annotations'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)


###################################
# 0) pre-set plotting aesthetics
in2mm<-25.4
my_theme = theme_classic(base_size = 6)

subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 'Astrocytes', 
               'Endothelial', 'Microglia', 'Mural', 'Oligos', 'Oligos_Pre', 
               'Interneuron')
othertypes_col = c(carto_pal(4, "Safe"), 
                   carto_pal(length(othertypes) -4 , "Vivid"))
names(othertypes_col) = othertypes
othertypes_col = othertypes_col[-c(1:4)]
othertypes_col = othertypes_col[c(7, 1:6)]
typecolors = c(subtypes_col, othertypes_col)

###################################
# 1) read in Logan snRNA dataset to plot
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'SCT') 

# set the groupings by the refined cell type column
## merge interneurons again
obj_merged$celltype3 = ifelse(grepl('Int', obj_merged$celltype3), 
                              'Interneuron',obj_merged$celltype3)
Idents(object = obj_merged) <- "celltype3"

# change the ordering of the cell types
levels(obj_merged) <- names(typecolors)
table(Idents(obj_merged))

####################################################################
# 2) plot the diagonal matrix of opioid receptors and their endogenous peptide

## plot the opioid receptors and ligands
markOPR1 = c('OPRM1', 'OPRD1', 'OPRK1', 'OPRL1') # opioid receptor genes
markOPR2 = c('PENK', 'PDYN', 'POMC', 'PNOC') # endorphin peptide genes
markerGenes = c(markOPR1,markOPR2 )

fig2_opioid_matrix_dotplot_fn = 
  here(PLOTDIR, 'plots', 'fig2_opioid_matrix_dotplot.all.pdf')

# https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
pdf(fig2_opioid_matrix_dotplot_fn, width = 80/in2mm, height =  70/in2mm)
DotPlot( obj_merged, features = markerGenes, cols = c("lightgrey", 'limegreen'),
  cluster.idents = F, scale = T, scale.by = "radius") +
  my_theme + scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  guides(size = FALSE) +
  theme(legend.position = 'bottom', plot.title= element_blank(),
      axis.text.x = element_text(angle = -30, vjust = 0, hjust=1),
      axis.title = element_blank(), 
      legend.key.size = unit(2, "mm"),
      legend.spacing.x = unit(1, 'mm'))
dev.off()



####################################################################
# 3) plot the diagonal matrix of opioid receptors and their endogenous peptide

## plot the opioid receptors and ligands
markerGenes2 = c('CLOCK', 'ARNTL', 'PER1', 'PER2', 'PER3', 'CRY1', 'CRY2',
                 'NFIL3', 'NPAS2', 'NR1D1', 'NR1D2', 'SIRT1', 'NAMPT')

fig2_circadian_matrix_dotplot_fn = 
  here(PLOTDIR, 'plots', 'fig2_circadian_matrix_dotplot.all.pdf')

# https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
pdf(fig2_circadian_matrix_dotplot_fn, width = 100/in2mm, height =  70/in2mm)
DotPlot( obj_merged, features = markerGenes2, cols = c("lightgrey", 'mediumblue'),
              cluster.idents = F, scale = T, scale.by = "radius") +
  my_theme + scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  theme(legend.position = 'bottom', plot.title= element_blank(),
        axis.text.x = element_text(angle = -30, vjust = 0, hjust=1),
        axis.title = element_blank(), 
        legend.key.size = unit(2, "mm"),
        legend.spacing.x = unit(1, 'mm'))
dev.off()

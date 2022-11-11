ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(ggsci)
library(dendextend) # install.packages('dendextend') # if don't have
library(ape)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## num parallel threads for Seurat
plan(sequential) # no parallel
options(future.globals.maxSize = 80e9)
library(here)

## set up the folders
DATADIR='data/raw_data'
PLOTDIR='figures/explanatory/figure2_the_celltype_annotations'

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
  LoadH5Seurat() 

# set the groupings by the refined cell type column
## merge interneurons again
obj_merged$celltype3 = ifelse(grepl('Int', obj_merged$celltype3), 
                              'Interneuron',obj_merged$celltype3)
obj_merged$celltype1 = ifelse(grepl('UNK_ALL', obj_merged$celltype1), 
                              obj_merged$celltype3, obj_merged$celltype1)
obj_merged$celltype1 = ifelse( obj_merged$celltype1 == 'Mural', 
                               'Mural/Fibroblast', obj_merged$celltype1)
table(obj_merged$celltype1)


# change the ordering of the cell types
Idents(object = obj_merged) <- "celltype1"
levels(obj_merged) <-c('MSNs', 'Interneurons', 'Astrocytes', 'Endothelial', 
                       'Microglia', 'Mural/Fibroblast', 'Oligos', 'Oligos_Pre')

## build the human cell type phylogenetic relationships
obj_merged <- BuildClusterTree(object = obj_merged)
tree_hg = Tool(object = obj_merged, slot = 'BuildClusterTree')
tree_hg = root(tree_hg, outgroup = 'MSNs')
class(tree_hg)


####################################
# 2) get the monkey cluster dendrogram of cell types
monkey_all = here('data/tidy_data/HeKleyman2021_macaque_striatum_data_processing', 
                  'rdas/GSE167920_Results_full_nuclei_processed_final.h5Seurat') %>% 
  LoadH5Seurat(assays = "integrated") %>% 
  subset(subset = region_name %in% c( "caudate", 'putamen')) %>% 
  RunPCA(verbose = FALSE) 

Idents(monkey_all) = 'cell_type_2'
levels(monkey_all) = levels(obj_merged)
table(Idents(monkey_all))

## make cluster
monkey_all = BuildClusterTree(monkey_all)
tree_rm = Tool(object = monkey_all, slot = 'BuildClusterTree')
tree_rm = root(tree_rm, outgroup = c('MSNs'))
class(tree_rm)

###################################
# 2) plot the diagonal matrix of reported marker genes & cell types
# Custom these kendo, and place them in a list
fig2_monkey_human_label_cocordance_tanglegram_fn = 
  here(PLOTDIR, 'plots', 'fig2_monkey_human_label_cocordance_tanglegram.all.pdf')

dl <- dendlist( tree_hg %>% as.dendrogram(), tree_rm %>% as.dendrogram())
pdf(fig2_monkey_human_label_cocordance_tanglegram_fn, width = 100/in2mm, height =  70/in2mm)
tanglegram(dl,  common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = TRUE, 
           highlight_branches_lwd=FALSE, margin_inner=7, lwd=2
)
dev.off()



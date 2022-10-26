## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## tree specific functions
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(tidytree)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/hierarchical_deg_analyses'
PLOTDIR='figures/exploratory/hierarchical_deg_analyses/plots'
dir.create(PLOTDIR, showWarnings = F, recursive = T)
dir.create(here(DATADIR, 'tables'), showWarnings = F, recursive = T)

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 8 cores
plan("multicore", workers = 8)
options(future.globals.maxSize = 80 * 1024^3)

##################################################
# 1) load in cell type labels for label transfer
## read in Logan BU snRNA dataset to label transfer
save_merged_fn = here('data/tidy_data/Seurat_projects', 
                      "OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat")

## load the integrated object for 
obj = save_merged_fn %>% LoadH5Seurat(assay = 'integrated')

cell_class = obj[[]] %>% as.data.frame() %>%
  filter(!duplicated(celltype3)) %>% 
  mutate(cell_class = case_when(grepl('NUD|Hyb|ICj', celltype3) ~ 'DNew',
                                grepl('^D1', celltype3) ~ 'D1',
                                grepl('^D2', celltype3) ~  'D2',
                                grepl('Int', celltype3) ~ 'INT',
                                TRUE ~ "GLIA")) %>%
  dplyr::select(celltype3, cell_class) %>% as.data.frame() %>% 
  split(x = .$celltype3,f = .$cell_class)

cell_class = cell_class[c('D1', 'D2', 'DNew', 'INT',"GLIA")]
cell_class_col = setNames(RColorBrewer::brewer.pal(5, 'Set1'), names(cell_class))

## grouping cell types by compartment
cell_class2 = obj[[]] %>% as.data.frame() %>%
  filter(!duplicated(celltype3)) %>% 
  mutate(cell_class = case_when(grepl('NUD|Hyb|ICj', celltype3) ~ 'DNew', 
                                grepl('Strio', celltype3) ~ 'Striosome', 
                                grepl('Matrix|Shell', celltype3) ~ 'Matrix',
                                grepl('Int', celltype3) ~ 'INT',
                                TRUE ~ "GLIA")) %>%
  dplyr::select(celltype3, cell_class) %>% as.data.frame() %>% 
  split(x = .$celltype3,f = .$cell_class)
cell_class2 = cell_class2[c('Striosome', 'Matrix','DNew',"INT", 'GLIA')]
cell_class2_col = setNames(RColorBrewer::brewer.pal(5, 'Set1'), names(cell_class2))

################################################################################
# 2) embedding in the dimensionality reduction, use integrated corrected version
dimRed  = Embeddings(obj, reduction = "pca")[,1:40] %>% as.data.frame()
head(dimRed)

## get the average high dim distance per cell type
dimRed = cbind(dimRed, obj[[]] %>% as.data.frame()) %>% 
  group_by(celltype3) %>%
  summarise_if(is.numeric, mean) %>% 
  column_to_rownames('celltype3') %>%
  select(starts_with('PC')) %>% as.matrix()


## make into a distance matrix and then a tree
Matrix <- as.matrix(t(dimRed))
distMat <- as.dist(1 - cor(Matrix, method = "pearson"))
set.seed(101)
my_tree <- ape::nj(distMat)
# my_tree = root(my_tree,outgroup =  'Microglia')
tree_dat = as.treedata(my_tree) %>% groupOTU(cell_class, "cell_class")

pdf(here(PLOTDIR, 'OUD_Striatum_refined_celltype.byMarker.pdf'), height = 4, width = 2)
ggtree(tree_dat, aes(color=cell_class)) + 
  geom_tiplab(align=TRUE, color = 'black') + 
  scale_color_manual(values = cell_class_col, name = 'Grouping') + 
  xlim(c(0,3.5)) + 
  guides(color = guide_legend(ncol = 1, title.position = 'top')) +
  theme(legend.position = 'bottom',
        legend.spacing.x = unit(.2, 'cm'), 
        legend.key.size = unit(.2, "cm"))
dev.off()

### make tree colored by compartment
Matrix <- as.matrix(t(dimRed))
distMat <- as.dist(1 - cor(Matrix, method = "pearson"))
set.seed(101)
my_tree <- ape::nj(distMat)
# my_tree = root(my_tree,outgroup =  'Microglia')
tree_dat2 = as.treedata(my_tree) %>% groupOTU(cell_class2, "cell_class")
pdf(here(PLOTDIR, 'OUD_Striatum_refined_celltype.byComp.pdf'), height = 4, width = 2)
ggtree(tree_dat2, aes(color=cell_class)) + 
  geom_tiplab(align=TRUE, color = 'black') + 
  scale_color_manual(values = cell_class2_col, name = 'Grouping') + 
  xlim(c(0,3.5)) + 
  guides(color = guide_legend(ncol = 1, title.position = 'top')) +
  theme(legend.position = 'bottom',
        legend.spacing.x = unit(.2, 'cm'), 
        legend.key.size = unit(.2, "cm"))
dev.off()

## export the data driven tree
tree_file = here(DATADIR, 'tables', 'OUD_Striatum_refined_celltype.byComp.nh')
phylo <- as.phylo(tree_dat2)
phylo$edge.length <- NULL
## print the newick text
write.tree(phylo, tree_file)
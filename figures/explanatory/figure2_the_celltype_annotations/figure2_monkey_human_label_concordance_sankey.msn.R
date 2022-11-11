ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(ggsci)
library(ggsankey) # devtools::instmsn_github("davidsjoberg/ggsankey")
library(here)


## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## num parmsnel threads for Seurat
plan(multicore, workers = 8) 
options(future.globals.maxSize = 80e9)

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
subtypes_col = subtypes_col[c(3,1,4,2, 5)]

###################################
# 1) read in Logan snRNA dataset to plot
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_msn_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat() 

# set the groupings by the refined cell type column
# change the ordering of the cell types
Idents(object = obj_merged) <- "celltype3"
levels(obj_merged) <- names(subtypes_col)


####################################
# 2) get the monkey cluster dendrogram of cell types
monkey_msn = here('data/tidy_data/HeKleyman2021_macaque_striatum_data_processing', 
                  'rdas/GSE167920_Results_MSNs_processed_final.h5Seurat') %>% 
  LoadH5Seurat(assays = "integrated") %>% 
  subset(subset = region_name %in% c( "caudate", 'putamen')) %>% 
  RunPCA(verbose = FALSE) 

## compute anchor cells b/t human snRNA dataset and Macaque dataset
anchors.mac_msn <- FindTransferAnchors(
  reference = monkey_msn, query = obj_merged, reduction = 'rpca',
  query.assay = 'integrated', reference.assay = 'integrated',
  normalization.method = 'SCT', dims = 1:30, reference.reduction = "pca")

## predict BU cell type w/ macaque 'cell_type_2' column
predictions.mac_msn <- TransferData(
  anchorset = anchors.mac_msn, refdata = monkey_msn$MSN_type, dims = 1:30, store.weights = F)
predictions.mac_msn = predictions.mac_msn %>% rownames_to_column('barcode') %>% 
  inner_join(obj_merged[[]] %>% rownames_to_column('barcode'))

## these are the cell-level monkey 
table(predictions.mac_msn$predicted.id)

predictions.mac_msn.long = predictions.mac_msn %>% 
  make_long(celltype3, predicted.id) %>% 
  mutate(node = factor(node, rev(names(subtypes_col))))

###################################
# 3) plot the diagonal matrix of human labels to monkey labels
fig2_monkey_human_label_cocordance_sankey_fn = 
  here(PLOTDIR, 'plots', 'fig2_monkey_human_label_cocordance_sankey.msn.pdf')

pdf(fig2_monkey_human_label_cocordance_sankey_fn, width = 40/in2mm, height =  67/in2mm)
ggplot(predictions.mac_msn.long, aes(x = x, next_x = next_x, node = node, 
                        next_node = next_node, fill = node, label = node)) +
  geom_sankey() + 
  geom_sankey_label(size = 1, color = "white", fill = "gray40") +
  scale_fill_manual(values = subtypes_col) +
  labs(x = NULL, y = NULL) + 
  scale_x_discrete(labels=c( "celltype1" = "Human\nLabels", 
                             "predicted.id" = "Rhesus\nMacaque\nLabels")) +
  theme_void(base_size = 6) + 
  theme(legend.position = 'none', axis.line=element_blank(),
        axis.text.y=element_blank(), axis.ticks.y=element_blank())  
dev.off()



library(AUCell)
library(readxl)
library(tidyverse)
library(here)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)
library(BiocParallel)
library(doParallel)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


DATADIR='data/tidy_data/AUCell_gene_set_activities'
dir.create(DATADIR, showWarnings = F)

######################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural/Fibroblast', 'Oligos', 'Oligos_Pre')
othertypes = c('Interneuron',  'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = c(carto_pal(length(othertypes) , "Vivid"))
names(othertypes_col) = othertypes

typecolors = c(subtypes_col, othertypes_col)

########################################
## 1) load in the refined cell type object
## read in Logan BU snRNA dataset to label transfer
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

obj_merged$celltype3 = ifelse(grepl('Int', obj_merged$celltype3), 'Interneuron',obj_merged$celltype3)
obj_merged$celltype3 = factor(obj_merged$celltype3, names(typecolors))
table(obj_merged$celltype3)
Idents(obj_merged) = 'celltype3'


meta = obj_merged[[]]
save_meta_fn = here(DATADIR, 'rdas', 'OUD_Striatum_refined_all_SeuratObj_N22.metadata.rds')
dir.create( here(DATADIR, 'rdas'), showWarnings = F)
saveRDS(meta, file=save_meta_fn)


########################################
## 2) read in th rhythmicity gene sets
alpha = 0.05
rhythmicity_fn = here('data/tidy_data/Xue2022_Rhythmicity_tables') %>% 
  list.files(pattern = 'obs_para', full.names = T)
names(rhythmicity_fn) = basename(rhythmicity_fn) %>% gsub('obs_para_|_biotype.csv', '', .)

geneSets = rhythmicity_fn %>% lapply(read.csv) %>% 
  lapply(function(x) x[x$pvalue < alpha, ] %>% pull('genes'))
lengths(geneSets)


## calculate the AUCell scores for each gene set
cells_AUC <- AUCell_run(obj_merged[["RNA"]]@data, geneSets, 
                        BPPARAM= MulticoreParam(20))

## save the AUCell scores and assignments
save_fn = here(DATADIR, 'rdas', 'AUCell_Xue_rhythmicity_scores_obs_para_refined_celltype_N22.rds')
dir.create( here(DATADIR, 'rdas'))
saveRDS(cells_AUC, file=save_fn)



######################################################
## 3) read in the OUD vs. CTL gain/loss of rhythmicity
rhythmicity_fn2 = here('data/tidy_data/Xue2022_Rhythmicity_tables') %>% 
  list.files(pattern = 'OUDvsCONT', full.names = T)
names(rhythmicity_fn2) = basename(rhythmicity_fn2) %>% gsub('_OUDvsCONT.csv', '', .)

geneSets2 = rhythmicity_fn2 %>% lapply(read.csv) %>% 
  lapply(function(x) 
    x %>% pivot_longer(cols = c('gainSig', 'lossSig'), 
                       names_to = 'group', values_to = 'isSig') %>% 
      filter(isSig))

geneSets2 = c(geneSets2 %>% lapply(function(x) x[x$group =='gainSig',] %>% pull(genes)), 
              geneSets2 %>% lapply(function(x) x[x$group =='lossSig',] %>% pull(genes)))
names(geneSets2) = paste(names(geneSets2), rep(c('gain', 'loss'), each = 2), sep = '_')
lengths(geneSets2)


## calculate the AUCell scores for each gene set
cells_AUC2 <- AUCell_run(obj_merged[["RNA"]]@data, geneSets2, 
                        BPPARAM= MulticoreParam(10))

## save the AUCell files
save_fn = here(DATADIR, 'rdas', 'AUCell_Xue_OUDvsCONT_gainLoss_rhythmicity_refined_celltype_N22.rds')
dir.create( here(DATADIR, 'rdas'))
saveRDS(cells_AUC2, file=save_fn)





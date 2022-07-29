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

## AUCell object for DNA damage estimates
library(AUCell)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

PLOTDIR='figures/exploratory/dna_damage'

#########################################
# 1) load in the DNA damage estimates
dnaDam_val = here('data/tidy_data/Seurat_projects/AUCell', 
                  '22_07_07_cells_AUC.rds') %>% readRDS() %>% 
  getAUC() %>% as.matrix() %>% t() %>% as.data.frame() 

dnaDam_class = here('data/tidy_data/Seurat_projects/AUCell', 
                    '22_07_07_cells_assignment.rds') %>% 
  readRDS() %>% unlist()


#####################################
# 2) load in cell type annotations 
obj = here('data/tidy_data/Seurat_projects', 
           "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

## ordered cell types
celltypes = obj[[]] %>% filter(!duplicated(celltype3)) %>% 
  mutate(class = case_when(grepl('^D[1-2]', celltype3) ~ 'MSN', 
                           grepl('^Int', celltype3) ~ 'INT', 
                           TRUE ~ 'GLIA'), 
         class = factor(class, c('MSN', 'INT', 'GLIA'))) %>% 
  arrange(class, celltype3) %>% pull(celltype3)

## load the unfiltered QC table
celltype_df = obj[[]] %>% dplyr::select(-contains('integrated_snn_res')) %>% 
  mutate(celltype3 = factor(celltype3, celltypes))

celltype_df = cbind(celltype_df, 'DNA_dam_val' = dnaDam_val[rownames(celltype_df),]) %>% 
  mutate(DNA_dam_class = ifelse(rownames(celltype_df) %in% dnaDam_class, 'damaged', 'undamaged')) %>% filter(!is.na(DNA_dam_val))

## check class is right
celltype_df %>% group_by(DNA_dam_class) %>% summarise(mean(DNA_dam_val, na.rm = T))




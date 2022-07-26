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

PLOTDIR='figures/exploratory/dna_damage'

##############################################################
# 1) load in cell type annotations 

## load the unfiltered QC table
celltype_df = here('data/tidy_data/tables',
             paste0("BU_OUD_Striatum_refined_all_SeuratObj_N22.txt.gz")) %>%
  read_tsv(show_col_types = FALSE) %>% dplyr::select(-contains('integrated_snn_res'))

celltypes = celltype_df %>% filter(!duplicated(celltype3)) %>% 
  mutate(class = case_when(grepl('^D[1-2]', celltype3) ~ 'MSN', 
                           grepl('^Int', celltype3) ~ 'INT', 
                           TRUE ~ 'GLIA'), 
         class = factor(class, c('MSN', 'INT', 'GLIA'))) %>% 
  arrange(class, celltype3) %>% pull(celltype3)


## load in the DNA damage estimates
dna_dam_df = here('data/tidy_data/Seurat_projects/AUCell',
                  '22_07_07_BU_OUD_Striatum_filtered_SCT_SeuratObj_N22_g_signature_assignments_master_auc_values.csv') %>% read_csv() %>% as.data.frame() %>% 
  full_join(celltype_df) %>% mutate(celltype3 = factor(celltype3, celltypes)) %>% 
  filter(!is.na(nCount_SCT))

head(dna_dam_df$DNA_dam_val)
tail(dna_dam_df$DNA_dam_val)
summary(dna_dam_df$DNA_dam_val)




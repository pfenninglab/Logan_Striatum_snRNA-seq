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
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F)
dir.create(here(PLOTDIR, 'rdas'), showWarnings = F)

###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH', 'Int-SST', 'Int-TH', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = carto_pal(length(othertypes), "Vivid")
names(othertypes_col) = othertypes


###########################################
# 1) load in the DNA damage estimates from 
# damaged vs. undamage - general DNA damage signature (for the glia)
# damaged vs. neurons - neuron specific DNA damage signature

## get the numerical DNA damage scores from 2 calculations
dnaDam_val_all = here('data/tidy_data/Seurat_projects/AUCell', 
                  '22_08_08_cells_AUC.rds') %>% readRDS() %>% 
  getAUC() %>% as.matrix() %>% t() %>% as.data.frame() %>% 
  dplyr::select(damaged_vs_undamaged)

dnaDam_val_neu = here('data/tidy_data/Seurat_projects/AUCell', 
                  '22_08_08_cells_AUC_neurons_only.rds') %>% readRDS() %>% 
  getAUC() %>% as.matrix() %>% t() %>% as.data.frame() %>% 
  dplyr::select(damaged_vs_undamaged_neurons)

## get the DNA damage classifications
dnaDam_class_all = here('data/tidy_data/Seurat_projects/AUCell', 
                    '22_08_08_cells_assignment.rds') %>% readRDS()
dnaDam_class_all_cutoff = dnaDam_class_all[['damaged_vs_undamaged']][[1]]$selected
dnaDam_class_all = dnaDam_class_all[['damaged_vs_undamaged']][[2]]
table(dnaDam_class_all %>% ss('_', 2))

dnaDam_class_neu= here('data/tidy_data/Seurat_projects/AUCell', 
                        '22_08_08_cells_assignment_neurons_only.rds') %>% readRDS()
dnaDam_class_neu_cutoff = dnaDam_class_neu[['damaged_vs_undamaged_neurons']][[1]]$selected
dnaDam_class_neu = dnaDam_class_neu[['damaged_vs_undamaged_neurons']][[2]]
table(dnaDam_class_neu %>% ss('_', 2))


#####################################
# 2) load in cell type annotations 
obj = here('data/tidy_data/Seurat_projects', 
           "OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
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

celltype_df = cbind(celltype_df, 
                    DNA_dam_val_all = dnaDam_val_all[rownames(celltype_df),], 
                    DNA_dam_val_neu = dnaDam_val_neu[rownames(celltype_df),]) %>% 
  mutate(DNA_dam_class_all = ifelse(rownames(celltype_df) %in% dnaDam_class_all, 
                                    'damaged', 'undamaged'),
         DNA_dam_class_neu = ifelse(rownames(celltype_df) %in% dnaDam_class_neu, 
                                    'damaged', 'undamaged'), 
         DNA_dam_class_neu = ifelse(is.na(DNA_dam_val_neu), NA, DNA_dam_class_neu)) %>% 
  filter(!is.na(DNA_dam_val_all) | !is.na(DNA_dam_val_neu))
  
table(celltype_df$DNA_dam_class_neu)

save_file = here(PLOTDIR, 'rdas', 'OUD_Striatum_refined_all_SeuratObj_N22.meta.DNAdam2.rds')
celltype_df %>% saveRDS(save_file)


## add to Seurat object meta.data
obj = obj[,rownames(celltype_df)]
obj$DNA_dam_val_all = celltype_df$DNA_dam_val_all
obj$DNA_dam_val_neu = celltype_df$DNA_dam_val_neu
obj$DNA_dam_class_all = celltype_df$DNA_dam_class_all
obj$DNA_dam_class_all = celltype_df$DNA_dam_class_all
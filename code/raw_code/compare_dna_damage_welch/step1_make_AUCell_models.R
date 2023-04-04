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
library(rcartocolor)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


DATADIR='data/tidy_data/compare_dna_damage_welch'
dir.create(DATADIR, showWarnings = F)
dir.create( here(DATADIR, 'rdas'))

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
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

obj_merged$celltype3 = ifelse(grepl('Int', obj_merged$celltype3), 'Interneuron',obj_merged$celltype3)
obj_merged$celltype3 = factor(obj_merged$celltype3, names(typecolors))
table(obj_merged$celltype3)
Idents(obj_merged) = 'celltype3'


########################################
## 2) read in th dna_dam gene sets
alpha = 0.05
stages = c('Stage1' = 'celltype_factor_NeuNplus_gh2ax_plus_vs_NeuNplus_gh2ax_minus_l2fc1_signature.csv', 
           'Stage2' = 'celltype_factor_NeuNminus_gh2ax_plus_vs_NeuNplus_gh2ax_minus_l2fc1_signature.csv')

dna_dam_fn = here('data/tidy_data/compare_dna_damage_welch/tables') %>% 
  list.files(pattern = 'celltype_factor', full.names = T)
names(dna_dam_fn) = names(stages)[match( basename(dna_dam_fn),stages)]

geneSets = dna_dam_fn %>% lapply(read.csv) %>% 
  lapply(function(x) x[x$padj < alpha, ] %>% pull('hg_symbol'))
lengths(geneSets)

## calculate the AUCell scores for each gene set in just the glia
obj_glia = obj_merged[,! obj_merged$celltype3 %in% 
                          c('D1-Matrix', 'D2-Matrix', 'D1-Striosome', 
                            'D2-Striosome','D1/D2-Hybrid', 'Interneuron')]

## calculate the AUCell scores for each gene set
glia_AUC <- AUCell_run(obj_glia[["RNA"]]@data, geneSets,BPPARAM= MulticoreParam(20))

## save the AUCell scores and assignments
save_fn = here(DATADIR, 'rdas', 'AUCell_Welch_dna_dam_stages_refined_gliatype_N22.rds')
saveRDS(glia_AUC, file=save_fn)



## calculate the AUCell scores for each gene set in just the neurons
obj_neuron = obj_merged[,obj_merged$celltype3 %in% 
                          c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 
                            'D2-Striosome','D1/D2-Hybrid', 'Interneuron')]

neuron_AUC <- AUCell_run(obj_neuron[["RNA"]]@data, geneSets, 
                        BPPARAM= MulticoreParam(20))

## save the AUCell scores and assignments
save_fn2 = here(DATADIR, 'rdas', 'AUCell_Welch_dna_dam_stages_refined_neurontype_N22.rds')
saveRDS(neuron_AUC, file=save_fn2)




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

DATADIR = 'data/tidy_data/AUCell_gene_set_activities'
here(DATADIR,c( 'loom', 'regulon')) %>% sapply(dir.create)

########################################
## 1) load in the filtered Seurat object
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

DefaultAssay(object = obj_merged) = 'RNA'
names(obj_merged[[]])
table(obj_merged$Case)
table(obj_merged$celltype3) %>% sort()

###################################################################
## 2) subset genes to just those analyzed in differential analyses
## this reduces compute by 33%
res.celltype = here('data/tidy_data/differential_expression_analysis/rdas', 
                    'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% readRDS()
genes = res.celltype[[1]] %>% pull(gene) %>% sort() %>% unique()
obj_subset = obj_merged[genes, ]

###################################################################
## 3) output the full dataset of cells w/ filtered genes to loom
out_loom_fn = here(DATADIR,'regulon', paste0("OUD_Striatum_refined_N22"))
out_loom <- as.loom(obj_subset, filename = out_loom_fn, verbose = T, overwrite = T)
out_loom$close_all() ## close 

###################################################################
## 4) downsample oligos to 2nd most populous cell type, astrocytes
## this reduces compute by 50%
ind_drop = which(obj_merged$celltype3 == 'Oligos') %>% sample(size = 54500)
obj_subset2 = obj_subset[, -ind_drop]
table(obj_subset2$Case)
table(obj_subset2$celltype3) %>% sort()

########################################
## 5) create loom object for each sample
for( sub in obj_merged$Case %>% sort() %>% unique()){
  out_loom_fn = here(DATADIR,'loom', paste0("OUD_Striatum_refined.",sub,".loom"))
  if(!file.exists(out_loom_fn)){
    print(paste('making loom for:', sub))
    obj_subset3 = obj_subset2 %>% subset(Case == sub)
    out_loom <- as.loom(obj_subset3, filename = out_loom_fn, verbose = T, overwrite = T)
    out_loom$close_all() ## close 
  }
}



#
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(tidyverse)
library(broom)
library(here)
library(rcartocolor)
library(lme4)
library(lmeTest)
library(tidymodels)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)
library(BiocParallel)
library(doParallel)

DATADIR='data/raw_data'
PROJDIR='figures/exploratory/explore_DRD3_genes'

## make data dirs
here(PROJDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

##########################################
## 1) load in the refined cell type object
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

## add the D3 normalized expression to the metadata
obj_merged$DRD3 = GetAssayData(object = obj_merged, assay = 'RNA', slot = "data")['DRD3',]

###################################################################
## 2) subset genes to just those analyzed in differential analyses
## this reduces compute by 33%
res.celltype = here('data/tidy_data/differential_expression_analysis/rdas', 
                    'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% readRDS()
genes = res.celltype[[1]] %>% pull(gene) %>% sort() %>% unique()
obj_subset = obj_merged[genes, ]


###################################################################
## 2) subset genes to just those analyzed in differential analyses
obj_subset$celltype3 %>% table()

obj_D1 = obj_subset[,obj_subset$celltype3 %in% c('D1-Matrix', 'D1-Striosome')]

mat_D1 = GetAssayData(object = obj_D1, assay = 'RNA', slot = "data")
inds = seq(nrow(mat_D1))
names(inds) = rownames(mat_D1)
cor_D1 = parallel::mclapply(inds, function(idx){
  cor.test( obj_D1$DRD3, mat_D1[idx,]) %>% tidy()
}, mc.cores = 10) %>% bind_rows(.id = 'gene')

cor_D1 = cor_D1 %>% mutate(
  FDR = p.adjust(p.value, 'fdr'), 
  p.bonf = p.adjust(p.value, 'bonferroni')) %>% 
  arrange(p.value)
sum(cor_D1$FDR < 0.05, na.rm = T)
sum(cor_D1$p.bonf < 0.05, na.rm = T)

out_d1 = here(PROJDIR,'tables', 'D1_MSN_genes_correlated_to_DRD3.xlsx')
writexl::write_xlsx(cor_D1, out_d1)



obj_D2 = obj_subset[,obj_subset$celltype3 %in% c('D2-Matrix', 'D2-Striosome')]
mat_D2 = GetAssayData(object = obj_D2, assay = 'RNA', slot = "data")
inds = seq(nrow(mat_D2))
names(inds) = rownames(mat_D2)

cor_D2 = parallel::mclapply(inds, function(idx){
  cor.test( obj_D2$DRD3, mat_D2[idx,]) %>% tidy()
}, mc.cores = 20) %>% bind_rows(.id = 'gene')

cor_D2 = cor_D2 %>% mutate(
  FDR = p.adjust(p.value, 'fdr'), 
  p.bonf = p.adjust(p.value, 'bonferroni')) %>% 
  arrange(p.value)
sum(cor_D2$FDR < 0.05, na.rm = T)
sum(cor_D2$p.bonf < 0.05, na.rm = T)

out_d2 = here(PROJDIR,'tables', 'D2_MSN_genes_correlated_to_DRD3.xlsx')
writexl::write_xlsx(cor_D2, out_d2)








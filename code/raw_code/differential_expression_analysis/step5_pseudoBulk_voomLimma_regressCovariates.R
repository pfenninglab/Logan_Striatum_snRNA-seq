## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(rlang)
library(writexl)

## main differential gene expression package
library(limma)
library(SingleCellExperiment)
library(edgeR)
library(swfdr)
library(jaffelab)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'

###################################
# 0) grab the phenotype data
df = here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds') %>% 
  readRDS() %>% colData() %>% as.data.frame() %>% 
  dplyr::select(ID:CaseIdx) %>% 
  mutate(celltype3 = celltype3 %>% make.names %>% factor(names(typecolors)),
         celltype3 = droplevels(celltype3), 
         celltype_class = case_when(grepl('^D|^Int', celltype3)~ 'Neuron',  
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia')))
df = df[colnames(z_clean), ]
indlist = split(rownames(df), droplevels(df$celltype3))

###################################
# 1) expression matrix to regress

v = readRDS(here(DATADIR, 'rdas', 'voomLimma_norm_object_N222pb.rds'))
designSV = readRDS(here(DATADIR, 'rdas', 'bigModelFitSVA_designMatrix.rds'))

## move the cell type by sdiagnoses before all the covariables
designSV = designSV %>% as.data.frame() %>% 
  dplyr::relocate(contains('celltype_dx'), .before = everything()) %>% as.matrix()
num_regress = sum(grepl('celltype_dx', colnames(designSV)))

## regress everything except the cell type by diagnoses & interaction effects
## the E slot is log2 CPM, https://rdrr.io/bioc/limma/man/EList.html
# y_clean = 2 ^ cleaningY(v$E, designSV, P = num_regress)
y_clean = cleaningY(v$E, designSV, P = num_regress)

## sanity checks on the regressed expression values
summary(y_clean['RBFOX3',]) # should be high
summary(y_clean['DRD1',]) # should be medium

saveRDS(y_clean, here(DATADIR, 'rdas', 'voomLimma_norm_object_N222pb_regressedCovariate.rds'))

## calculate gene-wise z-scores of samples, grouped by cell type
z_clean = lapply(indlist, function(i) {
  a = apply(y_clean[,i], 1, scale, center = TRUE, scale = TRUE)
  rownames(a) = i
  return(a)
}) %>% Reduce('rbind', .) %>% t()
z_clean = z_clean[,rownames(df)]
saveRDS(z_clean, here(DATADIR, 'rdas', 'voomLimma_zscore_object_N222pb_regressedCovariate.rds'))





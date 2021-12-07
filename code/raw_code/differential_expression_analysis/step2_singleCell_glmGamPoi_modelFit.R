## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(rlang)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## main differential gene expression package
library(glmGamPoi)
library(SingleCellExperiment)
library(DelayedArray)
library(HDF5Array)
library(Matrix.utils)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 8 cores
plan("multicore", workers = 8)
options(future.globals.maxSize = 20000 * 1024^2)

##################################################
# 1) load in cell type labels for label transfer
## read in Logan BU snRNA dataset to label transfer
h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
sce = loadHDF5SummarizedExperiment(h5Dir, prefix="BU_Run1_Striatum_merged_RNA_SeuratObj_N4")

###########################################################
# 2) Fit or get the glmGamPoi fit for cell type by DSM OUD dx
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_fit_fn = here(rdasDir, 'BU_Run1_Striatum_glmGamPoi_modelFit_N4.rds')

if(!file.exists(save_fit_fn)){
  ## fit the glmGamPoi model for single cell expression
  fit <- glm_gp(sce, design = ~ celltype2 + DSM.IV.OUD + DSM.IV.OUD:celltype2 - 1,
                size_factors = 'deconvolution', reference_level = "Oligos" , 
                on_disk = FALSE, verbose = TRUE)

  ## save the output of edgeRQLFDetRate differential state analyses
  saveRDS(fit, save_fit_fn)
} else{
  fit = readRDS(save_fit_fn)
}

#############################################################
# 3) test differential expression per cell type by DSM OUD dx
celltypes = c('D1-Matrix', 'D2-Matrix', 'D1/D2-Hybrid',  'D1-Striosome', 
              'D2-Striosome', 'Interneurons', 'Astrocytes', 'Endothelial',
              'Microglia', 'Mural/Fibroblast', 'Oligos', 'Oligos_Pre')
celltypes = setNames(celltypes, celltypes)

## for every cell type, compute the DEGs 
res <- lapply(celltypes, function(cell) 
  test_de(fit, contrast = `DSM.IV.OUDOUD`, subset_to = celltype2 == cell,
          full_design = ~ DSM.IV.OUD + log10(nCount_RNA) - 1,
          pseudobulk_by = orig.ident) )

####################################################################
## 4) save the output of edgeRQLFDetRate differential state analyses
res = lapply(res, function(x) x %>%  
               dplyr::rename('gene' = 'name', 'p_val' = 'pval', 'logFC' = 'lfc', 
                             'p_adj' = 'adj_pval', 'F' = 'f_statistic'))

rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'BU_Run1_Striatum_glmGamPoiRes_byCelltype2_N4.rds')
saveRDS(res, save_res_fn)






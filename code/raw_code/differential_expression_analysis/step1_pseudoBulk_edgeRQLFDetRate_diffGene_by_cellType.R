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
library(limma)
library(edgeR)

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
save_merged_fn = here('data/tidy_data/Seurat_projects', 
                        "BU_Run1_Striatum_filtered_SCT_SeuratObj_N4.h5Seurat")
## load only the scaled and raw RNA counts
obj = save_merged_fn %>% LoadH5Seurat(assay = 'RNA') 
head(obj[[]])

## convert to SingleCellExperiment for GlmGamPoi
sce = as.SingleCellExperiment(obj)

##############################################################
# 2) compute per-sample per cell type pseudobulk DGE profiles
## save to hdf5 file format for quick save/load
h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
saveHDF5SummarizedExperiment(sce, h5Dir, prefix="BU_Run1_Striatum_merged_RNA_SeuratObj_N4")
sce = loadHDF5SummarizedExperiment(h5Dir, prefix="BU_Run1_Striatum_merged_RNA_SeuratObj_N4")

## store cell type IDs (kids) and sample IDs (sids)
nk <- length(kids <- rlang::set_names(levels(factor(sce$celltype2))))
ns <- length(sids <- rlang::set_names(levels(factor(sce$orig.ident))))

## aggregate by cluster-sample to create pseudo-bulk count matrix
groups <- colData(sce)[, c("celltype2", "orig.ident")]
pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 

## split by cluster, transform & rename columns
pb <- split.data.frame(pb, rep(kids, ns)) %>% 
  lapply(function(u) magrittr::set_colnames(t(u), unname(sids)))

## create per-subject colData table
pd_colData = colData(sce) %>% as.data.frame()%>%
  dplyr::select(orig.ident, Region:Cause.of.Death, starts_with('DSM.IV')) %>%
  filter(!duplicated(orig.ident))
rownames(pd_colData) = pd_colData$orig.ident

## create SingleCellExperiment from pseudo bulk counts per cell type
(pb <- SingleCellExperiment(assays = pb, colData = pd_colData[colnames(pb[[1]]),]))

#############################################################################
# 3) use edgeRQLFDetRate to detect per-cell type differential state analyses 
## for ea. cell type, run edgeRQLFDetRate w/ default parameters
## https://www.nature.com/articles/nmeth.4612
## https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
res <- lapply(kids, function(k) {
  ## calculate the detection rate of genes in this cell type
  y <- assays(pb)[[k]]
  cdr <- scale(colMeans(y > 0))

  ## construct design & contrast matrix regressing out the DetRate
  design <- model.matrix(~ 0 + cdr.V1 + DSM.IV.OUD, data = cbind(colData(pb), cdr = cdr))
  contrast <- makeContrasts("DSM.IV.OUDOUD-DSM.IV.OUDCTL", levels = design)
  
  ## fit per-cell type edgeR QLR w/ DetRate
  y <- DGEList(y, remove.zeros = TRUE)
  y <- calcNormFactors(y)
  y <- estimateDisp(y, design = design)
  fit <- glmQLFit(y, design = design)
  fit <- glmQLFTest(fit, contrast = contrast)
  
  ## get top diff gene per cell type
  topTags(fit, n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = k) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR)  
})

####################################################################
## 4) save the output of edgeRQLFDetRate differential state analyses
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'BU_Run1_Striatum_edgeRQLFDetRateRes_byCelltype2_N4.rds')
saveRDS(res, save_res_fn)
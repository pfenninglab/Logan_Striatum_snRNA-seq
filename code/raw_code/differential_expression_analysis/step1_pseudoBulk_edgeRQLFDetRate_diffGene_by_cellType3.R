## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(rlang)
library(writexl)

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
## set to be parallel over 28 cores
plan("multicore", workers = 28)
options(future.globals.maxSize = 80 * 1024^3)

##################################################
# 1) load in cell type labels for label transfer
## read in Logan BU snRNA dataset to label transfer
save_merged_fn = here('data/tidy_data/Seurat_projects', 
                        "OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat")
## load only the scaled and raw RNA counts
obj = save_merged_fn %>% LoadH5Seurat(assay = 'RNA') 
head(obj[[]])

## convert to SingleCellExperiment for GlmGamPoi
sce = as.SingleCellExperiment(obj)

##############################################################
# 2) compute per-sample per cell type pseudobulk DGE profiles
## save to hdf5 file format for quick save/load
h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
saveHDF5SummarizedExperiment(sce, h5Dir, prefix="OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat", replace=TRUE)
sce = loadHDF5SummarizedExperiment(h5Dir, prefix="OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat")

## store cell type IDs (k.ids) and sample IDs (s.ids)
nk <- length(kids <- rlang::set_names(levels(factor(sce$celltype3))))
ns <- length(sids <- rlang::set_names(levels(factor(sce$orig.ident))))

## aggregate by cluster-sample to create pseudo-bulk count matrix
groups <- colData(sce)[, c("celltype3", "orig.ident")]
pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 

## split by cluster, transform & rename columns
split_names = gsub(paste(paste0('_', sids), collapse = '|'), '', rownames(pb) )
pb <- split.data.frame(pb, rlang::set_names(split_names)) %>% 
  lapply(function(u) magrittr::set_colnames(t(u), gsub(paste(paste0(kids, '_'), collapse = '|'), '', rownames(u) )))

head(names(pb))

## these cell types not found in every sample
sapply(pb, ncol)
exclude_celltypes = c('Mural', 'Int-CCK', 'Int-PTHLH', 'Int-SST')
pb = pb[!names(pb) %in% exclude_celltypes]

## create per-subject colData table
pd_colData = colData(sce) %>% as.data.frame()%>%
  dplyr::select(orig.ident, Region:Cause.of.Death, starts_with('DSM.IV')) %>%
  mutate(Pair = factor(Pair)) %>%
  filter(!duplicated(orig.ident))
rownames(pd_colData) = pd_colData$orig.ident

## create SingleCellExperiment from pseudo bulk counts per cell type
(pb <- SingleCellExperiment(assays = pb, colData = pd_colData[colnames(pb[[1]]),]))

#############################################################################
# 3) use edgeRQLFDetRate to detect per-cell type differential state analyses 
## for ea. cell type, run edgeRQLFDetRate w/ default parameters
## https://www.nature.com/articles/nmeth.4612
## https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
res <- lapply(kids[!kids %in% exclude_celltypes ], function(k) {
  ## calculate the detection rate of genes in this cell type
  y <- assays(pb)[[k]]
  cdr <- scale(colMeans(y > 0))

  ##check with ryan about analysis DEG analysis in Seney et al from George Tseng, R script
  ## construct design & contrast matrix regressing out the DetRate
  design <- model.matrix(~ 0 + cdr.V1 + Age + Sex + PMI + RIN + Region + DSM.IV.OUD, 
                         data = cbind(colData(pb), cdr = cdr))
  contrast <- makeContrasts("DSM.IV.OUDOUD", levels = design)
  
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
save_res_fn = here(rdasDir, 'OUD_Striatum_edgeRQLFDetRateRes_bycelltype3_N22.rds')
saveRDS(res, save_res_fn)

tablesDir =file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
save_res_fn2 = here(tablesDir, 'OUD_Striatum_edgeRQLFDetRateRes_bycelltype3_N22.xlsx')
names(res) = make.names(names(res))
res %>% lapply(function(x) x %>% arrange(p_adj)) %>% writexl::write_xlsx(save_res_fn2)


##############################
## 5) cursory look at the DEGs
sapply(res, function(x) sum(x$p_adj < 0.05, na.rm = T))
# Astrocytes    D1-Matrix D1-Striosome D1/D2-Hybrid    D2-Matrix D2-Striosome Interneurons    Microglia       Oligos   Oligos_Pre 
#         50          295           45           20          172           27            1           79          588           59 

lapply(res, function(x) x %>% arrange(p_adj) %>% head(10))
lapply(res, function(x) x[c('OPRK1', 'OPRM1', 'OPRD1', 'PDYN', 'PENK'),] %>% 
         filter(p_adj < 0.10) %>% arrange(p_adj) %>% head(10))

#lapply(res, function(x) x[c('Maddie Gene', 'OPRM1', 'OPRD1', 'PDYN', 'PENK'),] %>% 
         #filter(p_adj < 0.10) %>% arrange(p_adj) %>% head(10)
names(res)
res[['D1-Striosome']] %>% filter(p_adj < 0.05) %>% arrange(p_adj)

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
library(SingleCellExperiment)
library(DelayedArray)
library(HDF5Array)
library(Matrix.utils)
library(limma)
library(edgeR)

## to make the output data tidy
library(biobroom)

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
# 1) create or load pseudobulk sce object
save_pseudobulk =here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds')
if(!file.exists(save_pseudobulk)){
  ## load the single cell counts
  h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
  sce = loadHDF5SummarizedExperiment(h5Dir, prefix="BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat")
  
  ## merge interneurons again
  sce$celltype3 = ifelse(grepl('Int', sce$celltype3), 'Interneuron',sce$celltype3)
  table(sce$celltype3)
  
  ## aggregate by cluster-sample to create pseudo-bulk count matrix
  colData(sce)
  groups <- colData(sce)[, c("celltype3", "Case", 'Region')]
  pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
  dim(pb)
  
  ## split by cluster, transform & rename columns
  pb_colData = colData(sce) %>% as.data.frame() %>%
    rownames_to_column('match') %>% 
    mutate(Pair = factor(Pair), match = paste(celltype3, Case, Region, sep = '_')) %>% 
    filter(!duplicated(match)) %>% column_to_rownames('match')
  pb_colData = pb_colData[rownames(pb),]
  
  ## make sure this PD is correct
  with(pb_colData, table(Case, celltype3, Region))
  
  ## add number of cells per aggregate
  num_cells = groups %>% as.data.frame() %>% 
    mutate(tmp = paste(celltype3, Case, Region,sep= '_')) %>% 
    pull(tmp) %>% table()
  num_cells %>% as.numeric() %>% summary()
  pb_colData$numCells = num_cells[rownames(pb_colData)]
  
  ## add the gene detection rate
  pb_colData$cdr <- scale(rowMeans(pb > 0)) 
  
  ## create SingleCellExperiment from pseudo bulk counts across all cell types and region
  (pb <- SingleCellExperiment(assays = t(pb), colData = pb_colData))
  
  ## remap case index nested inside OUD dx
  remap_case_idx = split(pb$Case, pb$DSM.IV.OUD) %>% 
    lapply(function(x){
      x = setNames(LETTERS[as.numeric(factor(x))], x)
      x[!duplicated(x)]
    }) %>% unlist()
  names(remap_case_idx) = ss(names(remap_case_idx), '\\.', 2)
  pb$CaseIdx = remap_case_idx[as.character(pb$Case)]
  
  ## check this is correct
  table(pb$DSM.IV.OUD, pb$CaseIdx)
  table(pb$celltype3, pb$CaseIdx)
  table(pb$celltype3, pb$Case)
  
  saveRDS(pb, save_pseudobulk)
} else {
  pb = readRDS(save_pseudobulk)
}


####################################################
## 2) filter pseudobulk samples that have too few cells
pb = pb[, pb$numCells > 20]
pb$celltype3 = make.names(pb$celltype3)

####################################
# 3) normalization using voom-limma
y <- DGEList(counts = assays(pb)[[1]])

## filter out genes w/ low counts
A <- rowSums(y$counts)
isexpr <- A > 50
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y)

y <- calcNormFactors(y)

## construct design & contrast matrix regressing out the DetRate
design <- model.matrix(~ cdr + Age + Sex + PMI + RIN + Region + # co-variates 
                         celltype3 + # nested variability of cell type in subject
                         DSM.IV.OUD:celltype3 + # celltype-specific effect in OUD
                         DSM.IV.OUD:Region + # region-specific effect in OUD
                         DSM.IV.OUD:Sex + # sex-specific effect in OUD
                         DSM.IV.OUD, # main effect
                         data = colData(pb))

all.zero <- apply(design, 2, function(x) all(x==0))
all.zero[all.zero]
idx <- which(colnames(design) %in% c('celltype3Mural:CaseIdxE', 'celltype3Mural:DSM.IV.OUDOUD') | 
               all.zero)
design <- design[,-idx]

## voom precision weights and sample-wise quality weights normalization
v <- voomWithQualityWeights(y, design)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Case)
cor$consensus #0.2302863

## recalculate weights after adjusting for correlated samples from same subject
v <- voomWithQualityWeights(y, design, block = colData(pb)$Case, 
          correlation = cor$consensus)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Case)
cor$consensus #0.2296111

## fit the model
fit <- lmFit(v, design, block = colData(pb)$Case, correlation = cor$consensus)
fit <- eBayes(fit)

tests = decideTests(fit,method = "separate")
summary(tests)

save_fit = here(DATADIR, 'rdas', 'voomLimma_diffGene_bigModelFit.rds')
saveRDS(fit, file = save_fit)

save_voom = here(DATADIR, 'rdas', 'voomLimma_norm_object_N205pb.rds')
saveRDS(v, file = save_voom)

## get cell type DEGs
celltype_effect = tests %>% as.data.frame() %>% names() %>% 
  grep('celltype3', .,value = T) %>% grep('DSM.IV.OUDOUD', ., value = T) 
names(celltype_effect) = gsub('celltype3|:DSM.IV.OUDOUD', '', celltype_effect)
deg_celltypes = lapply(celltype_effect, topTable, fit = fit, n=Inf)

## other effects
main_effect = tests %>% as.data.frame() %>% names() %>% 
  grep('celltype3', .,value = T, invert = T) %>% grep('DSM.IV.OUDOUD', ., value = T) 
names(main_effect) = gsub('OUD$|:DSM.IV.OUDOUD', '', main_effect)
deg_main = lapply(main_effect, topTable, fit = fit, n=Inf)

deg_list = c(deg_main, deg_celltypes) %>% lapply(function(x){
  x %>% dplyr::rename('p_adj' = 'adj.P.Val') %>% rownames_to_column('gene')
})
sapply(deg_list, function(x) x[x$p_adj < 0.05,] %>% nrow())

####################################################################
## 4) save the output of voom_limma differential state analyses
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'BU_OUD_Striatum_voom_limmaRes_bigModel_N22.rds')
saveRDS(deg_list, save_res_fn)

tablesDir =file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
save_res_fn2 = here(tablesDir, 'BU_OUD_Striatum_voom_limmaRes_bigModel_N22.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(p_adj)) %>% writexl::write_xlsx(save_res_fn2)





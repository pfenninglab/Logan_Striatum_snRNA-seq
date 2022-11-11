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
save_pseudobulk =here(DATADIR, 'rdas', 'OUD_Striatum_refined_all_PseudoBulk_N22.sce.rds')
if(!file.exists(save_pseudobulk)){
  ## lo
  h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
  sce = loadHDF5SummarizedExperiment(h5Dir, prefix="OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat")
  
  ## merge interneurons again
  sce$celltype3 = ifelse(grepl('Int', sce$celltype3), 'Interneuron',sce$celltype3)
  table(sce$celltype3)
  
  ## store cell type IDs (k.ids) and sample IDs (s.ids)
  nk <- length(kids <- rlang::set_names(levels(factor(sce$celltype3))))
  ns <- length(sids <- rlang::set_names(levels(factor(sce$Case))))
  
  ## aggregate by cluster-sample to create pseudo-bulk count matrix
  colData(sce)
  groups <- colData(sce)[, c("celltype3", "Case")]
  pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
  dim(pb)
  
  ## split by cluster, transform & rename columns
  pb_colData = colData(sce) %>% as.data.frame() %>%
    rownames_to_column('match') %>% 
    mutate(Pair = factor(Pair), match = paste(celltype3, Case, sep = '_')) %>% 
    filter(!duplicated(match)) %>% column_to_rownames('match')
  pb_colData = pb_colData[rownames(pb),]
  
  ## make sure this PD is correct
  with(pb_colData, table(Case, celltype3))
  
  ## add number of cells per aggregate
  num_cells = groups %>% as.data.frame() %>% 
    mutate(tmp = paste(celltype3, Case,sep= '_')) %>% 
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

#############################################################################
# 3) use edgeRQLFDetRate to detect per-cell type differential state analyses 
## for ea. cell type, run edgeRQLFDetRate w/ default parameters
## https://www.nature.com/articles/nmeth.4612
## https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R

## construct design & contrast matrix regressing out the DetRate
design <- model.matrix(~ cdr + Age + Sex + PMI + RIN + # co-variates 
                         celltype3 + celltype3:CaseIdx  + # nested variability of cell type in subject
                         DSM.IV.OUD:celltype3 + # celltype-specific effect in OUD
                         DSM.IV.OUD:Sex + # sex-specific effect in OUD
                         DSM.IV.OUD, # main effect
                         data = colData(pb))

all.zero <- apply(design, 2, function(x) all(x==0))
all.zero[all.zero]
idx <- which(colnames(design) %in% c('celltype3Mural:CaseIdxE', 'celltype3Mural:DSM.IV.OUDOUD') | 
               all.zero)
design <- design[,-idx]


## fit per-cell type edgeR QLR w/ DetRate
y <- DGEList(assays(pb)[[1]], remove.zeros = TRUE)
y <- calcNormFactors(y)
y <- estimateDisp(y, design = design)
fit <- glmQLFit(y, design = design)

save_fit = here(DATADIR, 'rdas', 'edgeRQLFDetRate_diffGene_bigModelFit.rds')
saveRDS(fit, file = save_fit)

## get the big OUD vs. CTL effect
design2 = design
colnames(design2) = make.names(colnames(design2))

#########################################
## 4) get the main effect of OUD vs. CTL
contrast <- makeContrasts("DSM.IV.OUDOUD", levels = colnames(design2))
fit_main <- glmQLFTest(fit, contrast = contrast)

deg_main = topTags(fit_main, n = Inf, sort.by = "none")$table %>% 
  dplyr::mutate(gene = rownames(.), cluster_id = 'main') %>% 
  dplyr::rename(p_val = PValue, p_adj = FDR) %>% 
  arrange(p_val)

## get number of genes w/ OUD effect
alpha = 0.05
sum(deg_main$p_adj < alpha) # 1000


#########################################
## 5) get the sex effect of OUD vs. CTL
contrast <- makeContrasts("SexM.DSM.IV.OUDOUD", levels = colnames(design2))
fit_sex <- glmQLFTest(fit, contrast = contrast)

deg_sex = topTags(fit_sex, n = Inf, sort.by = "none")$table %>% 
  dplyr::mutate(gene = rownames(.), cluster_id = 'sex') %>% 
  dplyr::rename(p_val = PValue, p_adj = FDR) %>% 
  arrange(p_val)


## get number of genes w/ OUD effect interacting w/ Sex
sum(deg_sex$p_adj < alpha) # 892
head(deg_sex, 10)

#########################################
## 5) get the main effect of OUD vs. CTL
## get the contrasts & fits for all the cell types in the model
celltypes = unique(pb$celltype3) %>% sort() %>% rlang::set_names()
contrasts = sapply(celltypes, function(x) paste0("celltype3", x, ".DSM.IV.OUDOUD"))
celltypes = celltypes[contrasts %in% colnames(design2)]
contrasts = contrasts[contrasts %in% colnames(design2)] 
fit_list = lapply(contrasts, function(x) {
  contrast = makeContrasts(contrasts = x,levels = colnames(design2))
  glmQLFTest(fit, contrast = contrast)
  })

deg_celltype_list = lapply(celltypes, function(x){
  topTags(fit_list[[x]], n = Inf, sort.by = "none")$table %>% 
    dplyr::mutate(gene = rownames(.), cluster_id = x) %>% 
    dplyr::rename(p_val = PValue, p_adj = FDR) %>% 
    arrange(p_val)
})

sapply(deg_celltype_list, function(deg) sum(deg$p_adj < alpha))
# D1.D2.Hybrid    D1.Matrix D1.Striosome    D2.Matrix D2.Striosome  Endothelial 
#            1            5            2            7            1           16 
# Interneuron    Microglia       Oligos   Oligos_Pre 
#           1            5            9            2 
sapply(deg_celltype_list, function(deg) deg[deg$p_adj < alpha,'gene']) # DLGAP3 is interesting

#########################################
## 6) aggregate the DEGs 
## get the contrasts & fits for all the cell types in the model





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
library(sva)
library(swfdr)

## regress out the surrogate variables
library(jaffelab)

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
save_pseudobulk =here(DATADIR, 'rdas', 'OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds')
pb = readRDS(save_pseudobulk)


####################################################
## 2) filter pseudobulk samples that have too few cells
pb = pb[, pb$numCells > 20]
pb$celltype3 = make.names(pb$celltype3) %>% as.factor()
pb$Region = as.factor(pb$Region)
pb$Sex = as.factor(pb$Sex)


## construct design & contrast matrix regressing out the DetRate
design <- model.matrix(~ cdr + Age + Sex + PMI + RIN + Region + # co-variates 
                         celltype3 + # nested variability of cell type in subject
                         DSM.IV.OUD:celltype3 + # celltype-specific effect in OUD
                         DSM.IV.OUD:Sex + # sex-specific effect in OUD
                         DSM.IV.OUD:Region + # region-specific effect in OUD
                         DSM.IV.OUD, # main effect
                       data = colData(pb))

all.zero <- apply(design, 2, function(x) all(x==0))
all.zero[all.zero]
idx <- which(colnames(design) %in% c('celltype3Mural:DSM.IV.OUDOUD') | 
               all.zero)
design <- design[,-idx]


## construct the null model, used in regressing out the various factors
design0 <- model.matrix(~ cdr + Age + Sex + PMI + RIN + Region + celltype3, # co-variates 
                        data = colData(pb))


####################################
# 3) normalization using voom-limma
y <- DGEList(counts = assays(pb)[[1]])
dim(y)

## filter out genes w/ low counts
A <- rowMeans(y$counts)
isexpr <- A > 5
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y)

## filter out ribosomal genes, filter out mitochondria genes
keep.genes <- grep("^RP[SL]|^MT-",rownames(y), value = T, invert = T)
y = y[keep.genes, , keep.lib.size = FALSE]
dim(y) #  20273   205

# normalize counts
y <- calcNormFactors(y)

## voom precision weights and sample-wise quality weights normalization
v <- voomWithQualityWeights(y, design)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Case)
cor$consensus # 0.2732756

## recalculate weights after adjusting for correlated samples from same subject
v <- voomWithQualityWeights(y, design, block = colData(pb)$Case, 
          correlation = cor$consensus)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Case)
cor$consensus # 0.2719249


############################################################
# 4) Use surrogate variables to estimate unmodeled variation

## estimate the number of SVs from the adjusted
save_sva =here(DATADIR, 'rdas', 'OUD_Striatum_refined_all_PseudoBulk_N22.sva.rds')
if(! file.exists(save_sva)){
  (n.sv = num.sv(v$E, design, method="be", seed = set.seed(1))) #30
  svobj = sva(v$E, design, design0, n.sv=n.sv, B = 20)
  saveRDS(svobj, save_sva)
} else {
  svobj = readRDS(svobj)
}

## add the SVs to the model matrix
designSV = cbind(design, svobj$sv)
design0SV = cbind(design0, svobj$sv)

## recalculate sample quality weights after calculating the SVs
v <- voomWithQualityWeights(y, designSV, block = colData(pb)$Case, 
                            correlation = cor$consensus)
cor <- duplicateCorrelation(v, designSV, block = colData(pb)$Case)
cor$consensus # 0.08238471

save_voom = here(DATADIR, 'rdas', 'voomLimma_norm_object_N205pb.rds')
saveRDS(v, file = save_voom)

## regress out the un-modeled variation in the data
v_clean = v
v_clean$E <- cleaningY(v$E, designSV, P = ncol(design))

save_voom_clean = here(DATADIR, 'rdas', 'voomLimma_clean_object_N205pb.rds')
saveRDS(v_clean, file = save_voom_clean)

####################################################################
## 6) fit the model to get DEGs for various differences in the data
fit <- lmFit(v, designSV, block = colData(pb)$Case, correlation = cor$consensus)
fit <- eBayes(fit)

tests = decideTests(fit,method = "separate")
summary(tests)

save_fit = here(DATADIR, 'rdas', 'voomLimma_diffGene_bigModelFitSVA.rds')
saveRDS(fit, file = save_fit)

####################################################################
## 7) get the DEG tables and peak at interesting DEGs and no. of DEGs

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

## use SWFDR to increase power of detecting DEGs based on avg expression covariate
## https://pubmed.ncbi.nlm.nih.gov/30581661/
deg_list = c(deg_main, deg_celltypes) %>% lapply(function(x){
  x %>% dplyr::rename('p_adj' = 'adj.P.Val') %>% rownames_to_column('gene')
})

deg_list = deg_list %>% 
  lapply(function(x){
    x %>% mutate(p_adj =  lm_qvalue(P.Value, X=AveExpr)$q)
  })

####################################################################
## 8) save the output of voom_limma differential state analyses
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.rds')
saveRDS(deg_list, save_res_fn)

tablesDir =file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
save_res_fn2 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value)) %>% writexl::write_xlsx(save_res_fn2)

alpha = 0.01
# save tables DEGs w/ P.Value < alpha up and down regulated
save_res_fn3 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.lowConfCutOff.upReg.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value) %>% 
                      filter(P.Value < alpha, logFC > 0)) %>% 
  writexl::write_xlsx(save_res_fn3)

save_res_fn4 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.lowConfCutOff.dnReg.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value) %>% 
                      filter(P.Value < alpha, logFC < 0)) %>% 
  writexl::write_xlsx(save_res_fn4)


#######################################
## 9) check out some interesting genes

# FDR cutoff
sapply(deg_list, function(x) x[x$p_adj < 0.05,] %>% nrow())

# lower confidence cutoff
sapply(deg_list, function(x) x[x$P.Value < 0.01,] %>% nrow())


sapply(deg_list, function(x){
  x %>% filter(gene %in% c('CRHR1', 'LRP8'))%>% 
    pull(p_adj)
})

lapply(deg_list, function(x){
  x %>% filter(p_adj < alpha)%>% arrange(p_adj) %>% pull(gene) %>% 
    paste(collapse = ', ')
})

lapply(deg_list, function(x){
  x %>% filter(p_adj < alpha) %>% pull(gene) %>% 
    str_subset('^OPR|PDYN|PENK|POMC|^DRD') %>% 
    paste(collapse = ', ')
})


lapply(deg_list, function(x){
  x %>% filter(P.Value < 0.01) %>% pull(gene) %>% 
    str_subset('^OPR|PDYN|PENK|POMC|^DRD') %>% 
    paste(collapse = ', ')
})


lapply(deg_list, function(x){
  x %>% filter(p_adj < alpha) %>% pull(gene) %>% 
    str_subset('^FOX|^GAD|^CHAT$|^TH$|^KCN|^SCN|^SLC') %>% 
    paste(collapse = ', ')
})


lapply(deg_list, function(x){
  x %>% filter(P.Value < 0.01) %>% pull(gene) %>% 
    str_subset('^FOX|^GAD|^CHAT$|^TH$|^KCN|^SCN|^SLC') %>% 
    paste(collapse = ', ')
})

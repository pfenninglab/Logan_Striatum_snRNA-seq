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

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'

###################################################
## 1) read in the fitted big model w/ SVA features 
save_pseudobulk =here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds')
pb = readRDS(save_pseudobulk)
pb = pb[, pb$numCells > 20]
pb = pb[, pb$celltype3 != 'Mural'] # drop mural cells b/c too few

save_fit = here(DATADIR, 'rdas', 'voomLimma_diffGene_bigModelFitSVA.rds')
fit = readRDS(file = save_fit)

save_design = here(DATADIR, 'rdas', 'bigModelFitSVA_designMatrix.rds')
designSV = readRDS(file = save_design)

###########################################################
## 2) compute the differences b/t Dx within each cell type
celltypes=levels(factor(pb$celltype3 %>% make.names()))
designSV2 =designSV
colnames(designSV2) = make.names(colnames(designSV2))

## make the cell type contrasts
con_celltypes = sapply(setNames(celltypes, celltypes),function(cell) {
  OUD = paste0('celltype_dxOUD.',cell) %>% make.names()
  CTL = paste0('celltype_dxCTL.',cell) %>% make.names()
  
  OUDSexM = paste0('SexM.celltype_dxOUD.',cell) %>% make.names()
  CTLSexM = paste0('SexM.celltype_dxCTL.',cell) %>% make.names()

  paste('(',OUDSexM, '-',OUD, ')-(',CTLSexM, '-',CTL,')')
  })

## the first cell type x dx is the reference cell type
con_celltypes[1] = "( SexM.celltype_dxOUD.Astrocytes - celltype_dxOUD.Astrocytes )" 

## create the contrasts for OUD effect Between all cells or major classes
con_groups = c('All' = paste('(', paste(con_celltypes, collapse = ' + '), 
                             ')/', length(con_celltypes)), 
             'Neuron' = paste('(', paste(con_celltypes[grepl('^D|^Int',celltypes)], collapse = ' + '), 
                              ')/', sum(grepl('^D|^Int',celltypes))), 
             'Glia' =  paste('(', paste(con_celltypes[!grepl('^D|^Int',celltypes)], collapse = ' + '),
                             ')/', sum(!grepl('^D|^Int',celltypes))))

## refit the model based on these contrasts
cont.matrix <- makeContrasts(contrasts= c(con_groups, con_celltypes), levels=designSV2)
rownames(cont.matrix) = colnames(designSV)
fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()

## compute the DEGs from these contrasts
deg_list = lapply(setNames(colnames(cont.matrix),  names(c(con_groups, con_celltypes))), 
                  function(coef){
                    topTable(coef = coef, fit =fit2, n=Inf) %>% arrange(P.Value) %>% 
                      ## use SWFDR to increase power of detecting DEGs based on avg expression covariate
                      ## https://pubmed.ncbi.nlm.nih.gov/30581661/
                      mutate(adj.P.Val.Within =  lm_qvalue(P.Value, X=AveExpr)$q) %>%
                      dplyr::select(-adj.P.Val) %>% rownames_to_column('gene') 
                  })

## FDR correction Between all tests
deg_list = deg_list %>% data.table::rbindlist(idcol = 'celltype') %>% 
  mutate(adj.P.Val.Between =  lm_qvalue(P.Value, X=AveExpr)$q) %>%
  split(by = 'celltype')

# FDR cutoff
sapply(deg_list, function(x) x[x$adj.P.Val.Within < 0.05,] %>% nrow())
sapply(deg_list, function(x) x[x$adj.P.Val.Between < 0.05,] %>% nrow())

# lower confidence cutoff
sapply(deg_list, function(x) x[x$P.Value < 0.01,] %>% nrow())

####################################################################
## 3) save the output of voom_limma differential state analyses
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.sexInteraction.rds')
saveRDS(deg_list, save_res_fn)

tablesDir =file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
save_res_fn2 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.sexInteraction.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value)) %>% writexl::write_xlsx(save_res_fn2)

# save tables DEGs w/ P.Value < alpha up and down regulated
save_res_fn3 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.sexInteraction.lowConfCutOff.upReg.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value) %>% 
                      filter(P.Value < 0.01, logFC > 0)) %>% 
  writexl::write_xlsx(save_res_fn3)

save_res_fn4 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.sexInteraction.lowConfCutOff.dnReg.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value) %>% 
                      filter(P.Value < 0.01, logFC < 0)) %>% 
  writexl::write_xlsx(save_res_fn4)


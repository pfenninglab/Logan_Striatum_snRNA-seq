## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(rlang)
library(writexl)

## main differential gene expression package
library(SingleCellExperiment)
library(limma)
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
pb = pb[, pb$numCells > 15]
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

## make the cell type contrasts that will be the sex interaction
con_celltypes = sapply(setNames(celltypes, celltypes),function(cell) {
  cell = colnames(designSV2) %>% make.names() %>% 
    str_subset(paste0('\\.', cell, '\\.'))

  # get the female OUD and CTL groups
  OUD_F = cell %>% str_subset('OUD') %>% str_subset(paste0('\\.F\\.')) 
  CTL_F = cell %>% str_subset('CTL') %>% str_subset(paste0('\\.F\\.')) 
  N_CTL_F = length(CTL_F); N_OUD_F = length(OUD_F)
  
  OUD_F = OUD_F %>% paste(collapse = ' + ')
  CTL_F = CTL_F %>% paste(collapse = ' + ')
  
  # get the male OUD and CTL groups
  OUD_M = cell %>% str_subset('OUD') %>% str_subset(paste0('\\.M\\.')) 
  CTL_M = cell %>% str_subset('CTL') %>% str_subset(paste0('\\.M\\.'))
  N_OUD_M = length(OUD_M); N_CTL_M = length(CTL_M)

  OUD_M = OUD_M %>% paste(collapse = ' + ')
  CTL_M = CTL_M %>% paste(collapse = ' + ')
  
  # return the contrast, difference 
  paste('(',OUD_F,')/',N_OUD_F, '-(',CTL_F,')/',N_CTL_F, 
        '-(',OUD_M,')/',N_OUD_M, '+(',CTL_M,')/',N_CTL_M)
})

## proportion of each cell type
df_prop = read_tsv(here('data/tidy_data/tables/BU_OUD_Striatum_refined_celltype3_proportions.txt')) %>% 
  deframe()
names(df_prop) = names(df_prop) %>% make.names()
df_prop = df_prop[celltypes]

ind_neur = grepl('^D|^Int',celltypes)
ind_glia = !grepl('^D|^Int',celltypes)

## create the contrasts for OUD effect Between all cells or major classes
con_groups = c('All' = paste0('(', con_celltypes,')*', df_prop) %>% paste(collapse = ' + '), 
               'Neuron' = paste0('(', con_celltypes[ind_neur],')*', 
                                      df_prop[ind_neur]/sum(df_prop[ind_neur])) %>% paste(collapse = ' + '), 
               'Glia' =  paste0('(', con_celltypes[ind_glia],')*', 
                                     df_prop[ind_glia]/sum(df_prop[ind_glia])) %>% paste(collapse = ' + '))

## create the contrast matrix
con = c(con_groups, con_celltypes)
cont.matrix <- makeContrasts(contrasts= con, levels=designSV2)
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

####################################################################
## 8) save the output of voom_limma differential state analyses
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.SexInteraction.rds')
saveRDS(deg_list, save_res_fn)

tablesDir =file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
save_res_fn2 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.SexInteraction.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value)) %>% writexl::write_xlsx(save_res_fn2)

##############################################
## 9) merge w/ the sex-specific differences
save_res_fn2 = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds')
deg_df = deg_list %>% data.table::rbindlist() %>% 
  rename_at(vars(c(logFC:adj.P.Val.Between)), ~paste0(., '_Interaction'))

deg_list2 = readRDS(save_res_fn2) %>% data.table::rbindlist() %>% 
  inner_join(deg_df) %>% 
  relocate(
    ends_with('_SexF'), ends_with('_SexM'), starts_with('dir'), 
    ends_with('_Interaction'), .after = everything())  %>%
  split(by = 'celltype')
  

####################################################################
## 8) save the output of voom_limma differential state analyses
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.SexSpecific.rds')
saveRDS(deg_list2, save_res_fn)

tablesDir =file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
save_res_fn2 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.SexSpecific.xlsx')
deg_list2 %>% writexl::write_xlsx(save_res_fn2)


  

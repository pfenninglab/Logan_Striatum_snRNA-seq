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

## make the cell type contrasts
con_celltypes_male = sapply(setNames(celltypes, celltypes),function(cell) {
  cell = colnames(designSV2) %>% make.names() %>% 
    str_subset(paste0('\\.', cell, '\\.')) %>% str_subset(paste0('\\.M\\.'))
  OUD = cell %>% str_subset('OUD') 
  CTL = cell %>% str_subset('CTL')
  
  N_OUD = OUD %>% length()
  OUD = OUD %>% paste(collapse = ' + ')
  
  N_CTL= CTL %>% length()
  CTL = CTL %>% paste(collapse = ' + ')
  
  paste('(',OUD,')/',N_OUD, '-(',CTL,')/',N_CTL)
})

## proportion of each cell type
df_prop = read_tsv('data/tidy_data/tables/BU_OUD_Striatum_refined_celltype3_proportions.txt') %>% 
  deframe()
names(df_prop) = names(df_prop) %>% make.names()
df_prop = df_prop[celltypes]

ind_neur = grepl('^D|^Int',celltypes)
ind_glia = !grepl('^D|^Int',celltypes)

## create the contrasts for OUD effect Between all cells or major classes
con_groups_male = c('All' = paste0('(', con_celltypes_male,')*', df_prop) %>% paste(collapse = ' + '), 
               'Neuron' = paste0('(', con_celltypes_male[ind_neur],')*', 
                                 df_prop[ind_neur]/sum(df_prop[ind_neur])) %>% paste(collapse = ' + '), 
               'Glia' =  paste0('(', con_celltypes_male[ind_glia],')*', 
                                df_prop[ind_glia]/sum(df_prop[ind_glia])) %>% paste(collapse = ' + '))

## refit the model based on these contrasts
con_male = c(con_groups_male, con_celltypes_male)
con_female = con_male %>% str_replace_all('\\.M\\.', '.F.') %>% 
  setNames(names(con_male) %>% paste0('#SexF'))
names(con_male) = paste0(names(con_male), '#SexM')
cont.matrix <- makeContrasts(contrasts= c(con_female, con_male), levels=designSV2)
rownames(cont.matrix) = colnames(designSV)
fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()

## compute the DEGs from these contrasts
deg_list = lapply(setNames(colnames(cont.matrix), names(c(con_female, con_male))), 
                  function(coef){
                    topTable(coef = coef, fit =fit2, n=Inf) %>% arrange(P.Value) %>% 
                      ## use SWFDR to increase power of detecting DEGs based on avg expression covariate
                      ## https://pubmed.ncbi.nlm.nih.gov/30581661/
                      mutate(adj.P.Val.Within =  lm_qvalue(P.Value, X=AveExpr)$q) %>%
                      dplyr::select(-adj.P.Val) %>% rownames_to_column('gene') 
                  })

## FDR correction Between all tests
deg_list = deg_list %>% data.table::rbindlist(idcol = 'group') %>% 
  mutate(adj.P.Val.Between =  lm_qvalue(P.Value, X=AveExpr)$q, 
         celltype = ss(group, '#', 1), Sex = ss(group, '#', 2)) %>%
  dplyr::select(-group) %>% pivot_longer(-c(gene, celltype, Sex)) %>% 
  pivot_wider(id_cols = c(celltype, gene), names_from = c(name, Sex), 
              values_from = value) %>% 
  split(f = 'celltype')


####################################################################
## 3) save the output of voom_limma differential state analyses
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds')
saveRDS(deg_list, save_res_fn)

tablesDir = file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
save_res_fn2 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.xlsx')
deg_list %>% writexl::write_xlsx(save_res_fn2)

# # save tables DEGs w/ P.Value < alpha up and down regulated
# save_res_fn3 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.lowConfCutOff.upReg.xlsx')
# deg_list %>% lapply(function(x) x %>% arrange(P.Value) %>% filter(P.Value < 0.01, logFC > 0)) %>% 
#   writexl::write_xlsx(save_res_fn3)
# 
# save_res_fn4 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.lowConfCutOff.dnReg.xlsx')
# deg_list %>% lapply(function(x) x %>% arrange(P.Value) %>% filter(P.Value < 0.01, logFC < 0)) %>% 
#   writexl::write_xlsx(save_res_fn4)
# 

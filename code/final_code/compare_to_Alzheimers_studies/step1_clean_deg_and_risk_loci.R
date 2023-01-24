## packages for data table processing 
library(tidyverse)
library(data.table)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/compare_to_Alzheimers_studies'
dir.create(here(DATADIR, 'rdas'), recursive = T, showWarnings = F)

main_types = c('Neuron', 'Astrocytes', 'Microglia','Oligos',  'Oligos_Pre')
rename = c('Ast' = 'Astrocytes', 'Ex' = 'Neuron', 'In' = 'Neuron', 
           'Mic' = 'Microglia', 'Oli' = 'Oligos', 'Opc' = 'Oligos_Pre')

################################################################################
## 1) read in the Blanchard2022_41586_2022_5439_MOESM17_ESM.xlsx APOE-4 study DEGs

## Table S15 is snRNA_seq of APOE e4 vs. e3 alleles
df_list1 = here(DATADIR, 'tables', 'Blanchard2022_41586_2022_5439_MOESM17_ESM_S15.xlsx') %>% 
  readxl::read_xlsx() %>% 
  mutate(celltype = rename[celltype], comparison = 'APOE4_vs_APOE3_snRNA', 
         group = paste0(celltype, '#', comparison)) %>% 
  split(f = .$group) %>%  lapply(function(x){
    x %>%  arrange(P.Value) %>% filter(!duplicated(gene))
  })
names(df_list1)
saveRDS(df_list1, here(DATADIR, 'rdas', 'Blanchard2022_APOE_allele_celltype_DEG_list.rds'))

## Table S12 is RNA_seq of APOE e4 vs. e3 alleles in iPSC oligodendrocytes
df_list2= here(DATADIR, 'tables', 'Blanchard2022_41586_2022_5439_MOESM14_ESM_S12.csv') %>% 
  read_csv() %>% filter(status =='OK') %>% 
  mutate(celltype = 'Oligos', comparison = 'APOE4_vs_APOE3_iPSC', 
         group = paste0(celltype, '#', comparison)) %>% 
  split(f = .$group) %>%  lapply(function(x){
    x %>%  arrange(p_value) %>% filter(!duplicated(gene))
  })
names(df_list2)
saveRDS(df_list2, here(DATADIR, 'rdas', 'Blanchard2022_APOE_allele_iPSC_Oligo_DEG_list.rds'))

## Table S13 is RNA_seq of APOE e4 vs. e3 alleles in postmortem oligodendrocytes, by Wilcoxon
df_list3= here(DATADIR, 'tables', 'Blanchard2022_41586_2022_5439_MOESM15_ESM_S13.csv') %>% 
  read_csv() %>% dplyr::rename('gene' = 'feature') %>% 
  mutate(celltype = 'Oligos', comparison = make.names(grp), 
         comparison = gsub('\\.\\.\\.', '+', comparison),
         comparison = gsub('\\.\\.', '_', comparison),
         comparison = gsub('\\.only|\\.$', '', comparison),
         comparison = gsub('\\.vs\\.', '_vs_', comparison),
         comparison = gsub('APOE34\\+APOE44_vs_APOE33_AD\\.and\\.nonAD', 
                           'APOE4_vs_APOE3', comparison),
         comparison = paste0(comparison, '_PMB'),
         group = paste0(celltype, '#', comparison)) %>% 
  split(f = .$group) %>%  lapply(function(x){
    x %>% arrange(pval) %>% filter(!duplicated(gene))
  })
names(df_list3)
saveRDS(df_list3, here(DATADIR, 'rdas', 'Blanchard2022_APOE_PMB_Oligo_DEG_list.rds'))

##################################################################################
## 2) read in the Mathys2018_41586_2019_1195_MOESM4_ESM.xlsx AD pathology vs. non-path data
mathys_fn = here(DATADIR, 'tables', 'Mathys2018_41586_2019_1195_MOESM4_ESM.xlsx')
cols = c("gene", "IndModel.adj.pvals", "mean.1", 'mean.2', "IndModel.FC", 
         "MixedModel.z", "MixedModel.p", "DEGs.Ind.Model", "DEGs.Ind.Mix.models")
comparisons = c('ADPath_vs_NoPath', 'EarlyAD_vs_NoPath', 'LateAD_vs_EarlyAD')
cols2 = paste0(rep(cols, times = 3),'#' ,rep(comparisons, each = length(cols)))
     
## read in the supplemental table w/ the DEGs from Mathys et al per cell type          
df = lapply(set_names(names(rename)), readxl::read_excel, path = mathys_fn, 
                 skip = 2, col_names = cols2) %>% 
  rbindlist(idcol = 'celltype') %>% 
  mutate(celltype = rename[celltype]) 

## split up the 3 DEG tables by the different comparisons
df_list = list(df[,1:10], df[,c(1, 11:19)], df[,c(1, 20:28)])
names(df_list) = comparisons
df_list = df_list %>% lapply(function(x){ names(x) = c('celltype', cols); return(x)}) %>% 
  rbindlist(idcol = 'comparison') %>% filter(complete.cases(.)) %>% 
  mutate(group = paste0(celltype, '#', comparison)) %>% 
  dplyr::select(group, celltype, gene, mean.1, mean.2, IndModel.adj.pvals, IndModel.FC) %>% 
  split(f = .$group)

## merge the neuron types together w/ averaging
df_list[grepl('Neuron', names(df_list))] = df_list[grepl('Neuron', names(df_list))] %>% 
  lapply(function(x) x %>% arrange(MixedModel.p) %>% filter(!duplicated(gene)))
sapply(df_list, nrow)

saveRDS(df_list, here(DATADIR, 'rdas', 'Mathys2018_41586_2019_1195_MOESM4_ESM.rds'))


##################################################################################
## 3) read in the Wightman2021_41588_2021_921_MOESM4_ESM.xlsx GWAS of AD and ADBP
wightman_fn = here(DATADIR, 'tables', 'Wightman2021_41588_2021_921_MOESM4_ESM.xlsx')
df1 = readxl::read_excel(wightman_fn, sheet = 'Supplementary Table 10 Mapped G', 
                        skip = 1) %>% 
  dplyr::rename('gene' = 'Symbol') %>% 
  mutate(study = 'Wightman et al. 2021', case = 'AD & AD-by-proxy') %>% 
  dplyr::select(c(gene, study, case))

bellenguez_fn = here(DATADIR, 'tables', 'Bellenguez2022_41588_2022_1024_MOESM4_ESM.xlsx')
df2 = readxl::read_excel(bellenguez_fn, sheet = 'Supplementary Table 20', 
                         skip = 2) %>% 
  dplyr::rename('gene' = 'Gene') %>% 
  mutate(study = 'Bellenguez et al. 2022', case = 'AD & related dementias') %>% 
  dplyr::select(c(gene, study, case))

df = bind_rows(df1, df2)
saveRDS(df, here(DATADIR, 'rdas', 'AD_GWAS_genes_Wightman2021_Bellenguez2022.rds'))








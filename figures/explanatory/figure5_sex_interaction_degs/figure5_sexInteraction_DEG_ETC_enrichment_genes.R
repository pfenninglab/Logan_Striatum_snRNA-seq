ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggsci)
library(KEGGREST)
library(here)

DATADIR='data/tidy_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure5_sex_interaction_degs'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)


#################################################
# 1) load in the DEG table from the big analyses
alpha = 0.05
res_OUDwinSex = here(DATADIR,'differential_expression_analysis', 'rdas',
                     'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds') %>% 
  readRDS() %>% data.table::rbindlist() %>% 
  dplyr::filter(adj.P.Val.Between_SexM < alpha)

######################################
## 2) load in the genes 
pathways_df = msigdbr("human", category="C5") %>% as.data.table() %>% filter(gs_subcat == 'GO:CC')

keep = c( 'GO:0045283', 'GO:0042653', 'GO:0042652', 'GO:0045271',
          'GO:0045273', 'GO:0045275', 'GO:0045277', 'GO:0045281')
keep_match = 'fumarate|NADH|electron transport|cytochrome c|UBIQUINONE|respiratory chain|succinate|ATP synthase'

etc_pathways = pathways_df %>% 
  filter(grepl(keep_match, gs_description, ignore.case = T) | 
           grepl('RESPIRATORY_CHAIN_COMPLEX', gs_name)) %>% 
  arrange(gs_name) %>% split(x= .$human_gene_symbol, f = .$gs_name)
lengths(etc_pathways)

out = sapply(etc_pathways, function(x){
  sum(x %in% res_OUDwinSex$gene)/length(x)
}) 
out = out[out>0]
out/max(out)* 100 %>% round()

sapply(etc_pathways, function(x){
  x = res_OUDwinSex %>% arrange(gene) %>% filter(!duplicated(gene)) %>% 
    filter(gene %in% x) %>% pull(gene)
  paste(x, collapse = ',')
}) 

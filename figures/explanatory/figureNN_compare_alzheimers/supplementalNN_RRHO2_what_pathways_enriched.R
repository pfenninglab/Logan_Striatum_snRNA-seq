ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(data.table)
library(tidyverse)
library(RRHO2)
library(here)
library(msigdbr)

DATADIR='data/tidy_data/compare_to_Alzheimers_studies'

## make for this subdirs
PLOTDIR='figures/explanatory/figureNN_compare_alzheimers'
here(PLOTDIR, c('plots/RRHO2', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

###########################
# 1) load in the RRHO data
RRHO_list= readRDS(here(PLOTDIR,'rdas', paste0('compare_alzheimers_rrho.rds')))
search = c('ADPath_vs_NoPath', 'EarlyAD_vs_NoPath',
           'APOE34_vs_APOE33_nonAD_PMB', 'APOE4_vs_APOE3_PMB')
RRHO_list = RRHO_list[grepl(paste(search, collapse = '|'), names(RRHO_list))]

du_genes = lapply(RRHO_list, '[[', 'genelist_du') %>% lapply('[[', 'gene_list_overlap_du') ## the shared DEGs
ud_genes = lapply(RRHO_list, '[[', 'genelist_ud') %>% lapply('[[', 'gene_list_overlap_ud') ## the shared DEGs
lengths(c(ud_genes,du_genes))


#############################################
## 2) get gene ontologies, use human genes
## grab the H, Hallmark set of gene pathways
## grab the C2, which is the curated canonical pathway sets
## grab the C5, which is the Gene Ontology sets
## https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
pathways_df =  bind_rows(msigdbr("human", category="H"), 
                         msigdbr("human", category="C2"), 
                         msigdbr("human", category="C5"))

## get the SynGO gene ontologies, use human genes
syngo_df = readRDS(here('data/tidy_data/SynGO_bulk_download_release_20210225',
                        'rdas', 'syngo_annotations.rds'))
pathways_df = rbindlist(list(pathways_df, syngo_df), fill = T) 

## reshape/label for pathway naming purposes
pathways <-pathways_df %>% 
  mutate(
    gs_subcat = ifelse(is.na(gs_subcat) | gs_subcat == '', gs_cat, gs_subcat),
    gs_name = paste(gs_subcat, gs_name, sep ='#')) %>% 
  split(x = .$gene_symbol, f = .$gs_name)


## exclude the really really big gene sets
lengths(pathways) %>% summary()
pathways = pathways[lengths(pathways)<500]
length(pathways) # 21308

table(pathways_df$gs_cat)

pathways_df2 = pathways_df %>% 
  dplyr::select(gs_subcat, gs_name, gs_description) %>% distinct() %>% 
  dplyr::rename('pathway_group' = 'gs_subcat', 'pathway' = 'gs_name',
                'description' = 'gs_description')

###############################
## 3) overlap enrichment tests


ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(data.table)
library(tidyverse)
library(RRHO2)
library(here)

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



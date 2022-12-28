## conda activate r4
## packages for data table processing 
library(tidyverse)
library(data.table)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/SynGO_bulk_download_release_20210225'
dir.create(here(DATADIR, 'rdas'), recursive = T, showWarnings = F)

mappings = c("gs_name"='go_id',
             "gene_symbol"='hgnc_symbol',
             "gs_pmid"='pubmed_id',
             "gs_description"='go_name')

df = here(DATADIR,'tables', 'syngo_annotations.xlsx') %>% 
  readxl::read_xlsx() %>% 
  mutate(gs_cat='SynGO', 
         gs_url = 'https://www.syngoportal.org/data/SynGO_bulk_download_release_20210225.zip', 
         gs_exact_source = 'PMID:31171447') %>% 
  dplyr::rename(!!!mappings) %>% 
  dplyr::select(-c(id:hgnc_id,go_domain:evidence_nontracable_author_statement))

df %>% saveRDS(here(DATADIR,'rdas', 'syngo_annotations.rds'))
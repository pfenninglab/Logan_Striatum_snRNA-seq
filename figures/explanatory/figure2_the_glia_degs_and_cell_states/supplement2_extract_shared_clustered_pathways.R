## packages for data table processing 
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(data.table)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

FIGDIR='figures/explanatory/figure2_the_glia_degs_and_cell_states'
dir.create(here(FIGDIR, 'plots', 'gsea'), recursive = T, showWarnings = F)
dir.create(here(FIGDIR, 'tables'), recursive = T, showWarnings = F)

###################################
# 0) pre-set colors and cell types 
celltypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome',
              'D1/D2-Hybrid', 'Interneuron',  'Astrocytes', 'Endothelial', 
              'Microglia', 'Mural', 'Oligos', 'Oligos_Pre')

neuron_types = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome',
                 'D1/D2-Hybrid', 'Interneuron', 'Neuron') %>% make.names()

glia_types = c('Astrocytes', 'Endothelial', 'Microglia', 'Mural', 'Oligos', 
               'Oligos_Pre', 'Glia', 'All') %>% make.names()

in2mm<-25.4
alpha = 0.05

#################################################
# 1) load in the DEG table from the big analyses
save_res_fn = here('data/tidy_data/differential_expression_analysis', 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds')
res = readRDS(save_res_fn) %>% 
  lapply(function(x) x %>% filter(adj.P.Val.Between < alpha)) %>% 
  rbindlist() %>% filter(celltype %in% glia_types) %>% 
  arrange(adj.P.Val.Between)

# 2) read in the clustered pathway enrichment tables
curenrich_clustered = 
  here(FIGDIR, 'tables','figure2_clustered_glia_gsea_pathway_network.xlsx') %>% 
  readxl::read_xlsx() %>% mutate(
    leadingEdge = str_split(leadingEdge, ',')
  ) %>% unnest(leadingEdge) %>% dplyr::rename('gene'= 'leadingEdge') 

## interferon cluster
curenrich_clustered %>% filter(cluster_number == 3)  %>% 
  distinct(celltype, gene) %>% 
  inner_join(res) %>% arrange(gene)


## microglia cluster
curenrich_clustered %>% filter(celltype == 'Microglia')  %>% 
  distinct(celltype, gene) %>% 
  inner_join(res) %>% arrange(adj.P.Val.Between)

res %>% filter(celltype == 'Microglia') %>% filter(gene %in% c('ADGRB3'))
res %>% filter(gene %in% c('NFKBIA'))
res %>% filter(gene %in% c('APOE'))

## endothelial cluster
curenrich_clustered %>% filter(cluster_number %in% c(1, 2, 4, 5, 6)) %>%
  count(gene) %>% arrange(desc(n)) %>% filter(n>5) %>% 
  mutate(celltype = 'Endothelial') %>% inner_join(res)


## the down-regulated cluster
curenrich_clustered %>% 
  count(gene) %>% arrange(desc(n)) %>% filter(n>4) %>% inner_join(res) %>% 
  filter(!celltype %in% c('All', 'Glia')) %>% filter(grepl('GAB|GRM|GRIA|GRIN', gene))

curenrich_clustered %>% filter(cluster_number == 9)  %>% 
  count(gene) %>% arrange(desc(n)) %>% filter(n>=3) %>% pull(gene)

curenrich_clustered %>% filter(cluster_number == 9)  %>% 
  count(gene) %>% arrange(desc(n)) %>% filter(n>=3) %>% inner_join(res) %>% 
  filter(!celltype %in% c('All', 'Glia'))


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

FIGDIR='figures/explanatory/figure3_the_neuron_degs_and_cell_states'
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
  rbindlist() %>% filter(celltype %in% neuron_types) %>% 
  arrange(adj.P.Val.Between)

res %>% filter(gene == 'OGA')
res %>% filter(gene == 'MTRNR2L8')

# 2) read in the clustered pathway enrichment tables
curenrich_clustered = 
  here(FIGDIR, 'tables','figure3_clustered_gsea_pathway_network.xlsx') %>% 
  readxl::read_xlsx() %>% mutate(
    leadingEdge = str_split(leadingEdge, ',')
  ) %>% unnest(leadingEdge) %>% dplyr::rename('gene'= 'leadingEdge') 

## DNA replication & Cell Cycle checkpoint cluster
curenrich_clustered %>% filter(cluster_number %in% c( 2))  %>% 
  distinct(celltype, gene) %>% 
  inner_join(res) %>% arrange(gene)


## mitochondrial pathways
curenrich_clustered %>% filter(cluster_number %in% c(3, 4))  %>% 
  distinct(celltype, gene) %>% 
  inner_join(res) %>% arrange(gene)


## hypoxia pathways
curenrich_clustered %>% filter(cluster_number %in% c(11))  %>% 
  distinct(celltype, gene) %>% 
  inner_join(res) %>% arrange(gene)

## UV radiation pathways
curenrich_clustered %>% filter(cluster_number %in% c(10))  %>% 
  distinct(celltype, gene) %>% 
  inner_join(res) %>% arrange(celltype)


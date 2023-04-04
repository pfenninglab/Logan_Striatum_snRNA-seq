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

FIGDIR='figures/explanatory/figure5_sex_interaction_degs'
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
save_res_fn = here('data/tidy_data/differential_expression_analysis', 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds')
res = readRDS(save_res_fn) %>% rbindlist() %>% dplyr::select(-contains('dir')) %>% 
  pivot_longer(cols = -c(celltype, gene), names_sep = '_Sex', names_to = c('name', 'Sex')) %>% 
  pivot_wider(id_cols =c(celltype, gene, Sex)) %>% 
  filter(adj.P.Val.Between < alpha) %>% arrange(adj.P.Val.Between)

# 2) read in the clustered pathway enrichment tables
curenrich_clustered = 
  here(FIGDIR, 'tables','figure5_clustered_gsea_pathway_network_sexInteraction.xlsx') %>% 
  readxl::read_xlsx() %>% mutate(
    Sex = ss(OUD.v.CTL.in, 'Sex', 2),
    leadingEdge = str_split(leadingEdge, ',')
  ) %>% unnest(leadingEdge) %>% dplyr::rename('gene'= 'leadingEdge') 

## pull out FKBP5 gene
res %>% filter(gene == 'FKBP5')

## astrocyte pathways
curenrich_clustered %>% filter(cluster_number %in% c(6))  %>% 
  dplyr::filter(Sex == 'F') %>% distinct(celltype, gene, Sex) %>%
  inner_join(res) 

## men mitochondrial pathways
curenrich_clustered %>% filter(cluster_number %in% c(12, 15, 16))  %>% 
  dplyr::filter(Sex == 'M') %>% distinct(celltype, gene, Sex) %>%
  inner_join(res) %>% arrange(gene) %>% pull(gene)

curenrich_clustered %>% filter(cluster_number %in% c(15))  %>% 
  dplyr::filter(Sex == 'M') %>% distinct(celltype, gene, Sex) %>%
  inner_join(res) %>% arrange(gene)


## men slit-robo pathways
curenrich_clustered %>% filter(cluster_number %in% c(13))  %>% 
  dplyr::filter(Sex == 'M') %>% distinct(celltype, gene, Sex) %>%
  inner_join(res) %>% arrange(gene)

## UV radiation pathways
curenrich_clustered %>% filter(cluster_number %in% c(10))  %>% 
  distinct(celltype, gene) %>% 
  inner_join(res) %>% arrange(gene)
## packages for data table processing 
library(tidyverse)
library(igraph)
library(network)
library(sna)
library(stringr)
library(RColorBrewer)
library(data.table)
library(here)
library(leiden)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
source(here('code/final_code/Rutils/igraph_pathway_clustering.R'))

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

###########################################
# 1) read in the pathway enrichment tables
xl_path = here('data/tidy_data/differential_expression_analysis', 
               'tables/OUD_Striatum_GSEA_enrichment_msigdb_H_C2_C5_SynGO.xlsx')
readxl::excel_sheets(path = xl_path)

## read in the comparison of OUD vs. CTL averaged across all subtypes
gsea_df <- readxl::read_excel(path = xl_path, sheet = 'All') %>% 
  filter(!celltype %in% c('All', 'Glia', 'Neuron')) %>% 
  mutate(celltype = factor(celltype, c(neuron_types, glia_types)))


pdf(here(FIGDIR, 'plots', 'gsea', 'sfig_nes_by_pval.pdf'), 
    height = 5, width = 7)
ggplot(gsea_df, aes(x = abs(NES), y = -log10(pval) )) + 
  geom_smooth(method = 'lm', se = F, color = 'black', alpha = 0.5) +
  geom_point(aes(size = size), alpha = .3, pch = 20) +
  facet_wrap(~celltype) + theme_bw(base_size = 10)
dev.off()

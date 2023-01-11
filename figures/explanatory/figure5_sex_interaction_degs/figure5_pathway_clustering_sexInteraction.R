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

"%ni%" <- function(x, table) !(match(x, table, nomatch = 0) > 0)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)
source(here('code/final_code/Rutils/igraph_pathway_clustering.R'))

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

###########################################
# 1) read in the pathway enrichment tables
xl_path = here('data/tidy_data/differential_expression_analysis', 
               'tables/OUD_Striatum_GSEA_enrichment_msigdb_H_C2_C5_SynGO.xlsx')
readxl::excel_sheets(path = xl_path)

## read in the comparison of OUD vs. CTL averaged across all subtypes
interaction_df = readxl::read_excel(path = xl_path, sheet = 'SexInteraction')
gsea_female = readxl::read_excel(path = xl_path, sheet = 'SexF')
gsea_male = readxl::read_excel(path = xl_path, sheet = 'SexM')

## filter out the pathways are in the same direction
gsea_female2 = gsea_female %>% 
  mutate(NES_other = gsea_male$NES[match(pathway, gsea_male$pathway)]) %>% 
  filter(pathway %ni% gsea_male$pathway | sign(NES_other) != sign(NES)) %>% 
  dplyr::select(-NES_other)

gsea_male2 = gsea_male %>% 
  mutate(NES_other = gsea_female$NES[match(pathway, gsea_female$pathway)]) %>% 
  filter(pathway %ni% gsea_female$pathway | sign(NES_other) != sign(NES)) %>% 
  dplyr::select(-NES_other)

## subset by female and male NES to create two pathway networks
with(gsea_female2, table(celltype, OUD.v.CTL.in))
net_female=make_igraph_from_pathways(gsea_female2)
V(net_female)$shape = ifelse(V(net_female)$pathway %in% interaction_df$pathway, 'square', 'circle')

with(gsea_male2, table(celltype, OUD.v.CTL.in))
net_male=make_igraph_from_pathways(gsea_male2)
V(net_male)$shape = ifelse(V(net_male)$pathway %in% interaction_df$pathway, 'square', 'circle')


##################################################
# 2)  add the colors to the NES of either pathways
nes_max = max(c(V(net_female)$NES, V(net_male)$NES)) 
nes_min = min(c(V(net_female)$NES, V(net_male)$NES))
nes_range = nes_max - nes_min
  
col = circlize::colorRamp2(seq(nes_min, nes_max, nes_range/10), brewer.pal(11,"PiYG"))
V(net_female)$color= col(V(net_female)$NES)
V(net_male)$color= col(V(net_male)$NES)


##############################################################################
## trim and cluster the female network; CHANGE target_degree if you have too 
## many edges, clustering will not work + you will get a hairball network
net_female.sp = trim_edges(net_female, target_degree = .05)
clp_female <- cluster_leiden(net_female.sp, resolution_parameter = 1.5, 
                      objective_function = 'modularity', n_iterations = 10)
pdf(here(FIGDIR, 'plots','gsea','figure5_clustered_female_gsea_pathway_network.pdf'))
plot(clp_female, net_female.sp, vertex.label=NA)
dev.off()

net_female_list = trim_clustered_nodes(net_female.sp, clp_female, min_nodes = 5)
V(net_female_list$net)$hub_score = hub_score(net_female_list$net)$vector
V(net_female_list$net)$label.cex = .5

##############################################################################
## trim and cluster the male network; CHANGE target_degree if you have too 
## many edges, clustering will not work + you will get a hairball network
net_male.sp = trim_edges(net_male, target_degree = .05)
clp_male <- cluster_leiden(net_male.sp, resolution_parameter = .08, 
                      objective_function = 'modularity',  n_iterations = 10)
pdf(here(FIGDIR, 'plots','gsea','figure5_clustered_male_gsea_pathway_network.pdf'))
plot(clp_male, net_male.sp, vertex.label=NA)
dev.off()

net_male_list = trim_clustered_nodes(net_male.sp, clp_male, min_nodes = 4)
V(net_male_list$net)$mem = V(net_male_list$net)$mem + max(V(net_female_list$net)$mem)
V(net_male_list$net)$hub_score = hub_score(net_male_list$net)$vector
V(net_male_list$net)$label.cex = .5


###################################################################
## clean up the labels, merging the pathways from both clusterings
vert.meta =  bind_rows(vertex_attr(net_female_list$net), 
                       vertex_attr(net_male_list$net)) %>% as.data.frame()
to_label = vert.meta %>% clean_pathways()
to_label_num = vert.meta %>% dplyr::select(name, mem) %>% deframe()

## make sets of 3 plots for the current female network
prefix_female = here(FIGDIR, 'plots','gsea','figure5_clustered_female_gsea_pathway_network')
plot_clusterings_3set(net_female_list, prefix_female, to_label, to_label_num, 
                      height = 50/in2mm, width = 50/in2mm)
plot_pt_legend(V(net_male_list$net)$size, prefix_female, height = 20/in2mm, width = 10/in2mm)
plot_col_legend(prefix = prefix_female, col = col, fontsize = 5,
                height = 8/in2mm, width = 30/in2mm)


## make sets of 3 plots for the current male network
prefix_male = here(FIGDIR, 'plots','gsea','figure5_clustered_male_gsea_pathway_network')
plot_clusterings_3set(net_male_list, prefix_male, to_label, to_label_num, 
                      height = 50/in2mm, width = 50/in2mm)


## add the clustered pathways back into the 
to_label_num_female = vertex_attr(net_female_list$net) %>% as.data.frame() %>% 
  dplyr::select(name, mem) %>% deframe()
to_label_num_male = vertex_attr(net_male_list$net) %>% as.data.frame() %>% 
  dplyr::select(name, mem) %>% deframe()

curenrich_clustered = bind_rows(
  gsea_female2 %>% mutate(cluster_number = to_label_num_female[ paste0(celltype, '#', pathway)]),
  gsea_male2 %>% mutate(cluster_number = to_label_num_male[ paste0(celltype, '#', pathway)] ),
  interaction_df %>% mutate(cluster_number = to_label_num[ paste0(celltype, '#', pathway)]))  %>% 
  arrange(cluster_number, padj) %>% distinct() %>% 
  relocate(cluster_number, pathway, description,.after = 'celltype') %>% 
  writexl::write_xlsx(here(FIGDIR, 'tables','figure5_clustered_gsea_pathway_network_sexInteraction.xlsx'))

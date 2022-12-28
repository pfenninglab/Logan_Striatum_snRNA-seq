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

FIGDIR='figures/explanatory/figure3_the_degs_and_cell_states'
dir.create(here(FIGDIR, 'plots', 'gsea'), recursive = T, showWarnings = F)
dir.create(here(FIGDIR, 'tables'), recursive = T, showWarnings = F)


###################################
# 0) pre-set colors and cell types 
celltypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome',
              'D1/D2-Hybrid', 'Interneuron',  'Astrocytes', 'Endothelial', 
              'Microglia', 'Mural', 'Oligos', 'Oligos_Pre')

neuron_types = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome',
                 'D1/D2-Hybrid', 'Interneuron', 'Neuron') %>% make.names()
type_shape = setNames(names(igraph:::.igraph.shapes)[-7][seq_along(neuron_types)],
                      neuron_types)
in2mm<-25.4

###########################################
# 1) read in the pathway enrichment tables
xl_path = here('data/tidy_data/differential_expression_analysis', 
               'tables/OUD_Striatum_GSEA_enrichment_msigdb_H_C2_C5_SynGO.xlsx')
readxl::excel_sheets(path = xl_path)

## read in the comparison of OUD vs. CTL averaged across all subtypes
gsea_list <- readxl::read_excel(path = xl_path, sheet = 'All')
gsea_df = gsea_list%>% filter(celltype %in% neuron_types)
  
## subset by positive and negative NES to create two pathway networks
with(gsea_df, table(celltype, NES > 0 ))
net_pos=make_igraph_from_pathways(gsea_df %>% filter(NES > 0))
net_neg=make_igraph_from_pathways(gsea_df %>% filter(NES < 0))

#########################################################
# 2) craft the networks to create the pathway clustering

## add the colors to the POSITIVE NES
pos_range = max(V(net_pos)$NES) - min(V(net_pos)$NES)
col_pos = circlize::colorRamp2(seq(min(V(net_pos)$NES), max(V(net_pos)$NES), pos_range/8), 
                               brewer.pal(9,"BuGn"))
V(net_pos)$color= col_pos(V(net_pos)$NES)

## add the colors to the NEGATIVE NES
neg_range = max(V(net_neg)$NES) - min(V(net_neg)$NES)
col_neg = circlize::colorRamp2(seq(min(V(net_neg)$NES), max(V(net_neg)$NES), neg_range/8), 
                               rev(brewer.pal(9,"PuRd")))
V(net_neg)$color= col_neg(V(net_neg)$NES)

##############################################################################
## trim and cluster the POSITIVE network; CHANGE target_degree if you have too 
## many edges, clustering will not work + you will get a hairball network
net_pos.sp = trim_edges(net_pos, target_degree = .15)
clp_pos <- cluster_leiden(net_pos.sp, resolution_parameter = 1, 
                      objective_function = 'modularity', n_iterations = 10)
pdf(here(FIGDIR, 'plots','gsea','figure3_clustered_pos_gsea_pathway_network.pdf'))
plot(clp_pos, net_pos.sp, vertex.label=NA)
dev.off()

net_pos_list = trim_clustered_nodes(net_pos.sp, clp_pos, min_nodes = 5)
V(net_pos_list$net)$hub_score = hub_score(net_pos_list$net)$vector
V(net_pos_list$net)$label.cex = .5

##############################################################################
## trim and cluster the NEGATIVE network; CHANGE target_degree if you have too 
## many edges, clustering will not work + you will get a hairball network
net_neg.sp = trim_edges(net_neg, target_degree = .25)
clp_neg <- cluster_leiden(net_neg.sp, resolution_parameter = .5, 
                      objective_function = 'modularity',  n_iterations = 10)
pdf(here(FIGDIR, 'plots','gsea','figure3_clustered_neg_gsea_pathway_network.pdf'))
plot(clp_neg, net_neg.sp, vertex.label=NA)
dev.off()

net_neg_list = trim_clustered_nodes(net_neg.sp, clp_neg, min_nodes = 3)
V(net_neg_list$net)$mem = V(net_neg_list$net)$mem + max(V(net_pos_list$net)$mem)
V(net_neg_list$net)$hub_score = hub_score(net_neg_list$net)$vector
V(net_neg_list$net)$label.cex = .5


###################################################################
## clean up the labels, merging the pathways from both clusterings
vert.meta =  bind_rows(vertex_attr(net_pos_list$net), 
                       vertex_attr(net_neg_list$net)) %>% as.data.frame()
to_label = vert.meta %>% clean_pathways()
to_label_num = vert.meta %>% dplyr::select(name, mem) %>% deframe()

## make sets of 3 plots for the current POSITIVE network
prefix_pos = here(FIGDIR, 'plots','gsea','figure3_clustered_pos_gsea_pathway_network')
plot_clusterings_3set(net_pos_list, prefix_pos, to_label, to_label_num, 
                      height = 60/in2mm, width = 60/in2mm)
plot_pt_legend(V(net_neg_list$net)$size, prefix_pos, height = 20/in2mm, width = 10/in2mm)
plot_col_legend(prefix = prefix_pos, col = col_pos, fontsize = 5,
                height = 8/in2mm, width = 30/in2mm)


## make sets of 3 plots for the current NEGATIVE network
prefix_neg = here(FIGDIR, 'plots','gsea','figure3_clustered_neg_gsea_pathway_network')
plot_clusterings_3set(net_neg_list, prefix_neg, to_label, to_label_num, 
                      height = 60/in2mm, width = 60/in2mm)
plot_col_legend(prefix = prefix_neg, col = col_neg, fontsize = 5,
                height = 8/in2mm, width = 30/in2mm)


## add the clustered pathways back into the 
curenrich_clustered = gsea_df %>% 
  mutate(cluster_number = to_label_num[ paste0(celltype, '#', pathway)]) %>% 
  arrange(cluster_number, padj) %>%
  relocate(cluster_number, pathway, description,.after = 'celltype') %>% 
  writexl::write_xlsx(here(FIGDIR, 'tables','figure3_clustered_gsea_pathway_network.xlsx'))

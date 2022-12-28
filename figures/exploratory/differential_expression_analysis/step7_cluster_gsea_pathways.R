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

FIGDIR='figures/exploratory/differential_expression_analysis'
dir.create(here(FIGDIR, 'plots'), recursive = T, showWarnings = F)
dir.create(here(FIGDIR, 'tables'), recursive = T, showWarnings = F)


###################################
# 0) pre-set colors and cell types 
celltypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome',
                'D2-Striosome','D1/D2-Hybrid', 
              'Interneuron',  'Astrocytes', 'Endothelial', 
              'Microglia', 'Mural', 'Oligos', 'Oligos_Pre')

#################################################################
#example of getting nodes and edges from RERconverge enrichment:
xl_path = here('data/tidy_data/differential_expression_analysis', 
               'tables/OUD_Striatum_GSEA_enrichment_msigdb_H_C2_C5.xlsx')
tab_names <- readxl::excel_sheets(path = xl_path)

gsea_list <- lapply(set_names(tab_names), function(x) readxl::read_excel(path = xl_path, sheet = x))
gsea_df = gsea_list %>% rbindlist() 
  
with(gsea_df, table(OUD.v.CTL.in,celltype, NES > 0 ))
with(gsea_df, table(OUD.v.CTL.in,celltype))

## subset to just the main comparisons
curenrich = gsea_df %>% filter(NES>0) %>% 
  filter(OUD.v.CTL.in == 'All' & celltype != 'All') %>% 
  # filter(grepl('^D|In|Neu', celltype)) %>% 
  arrange(padj) %>% filter(!duplicated(pathway)) 
net=make_igraph_from_pathways(curenrich)
type_shape = c('Neuron' = "vrectangle", 'Glia' = 'circle')

##################################
#create color mapping from node values:
val_range = max(V(net)$NES) - min(V(net)$NES)
mycol = circlize::colorRamp2(seq(min(V(net)$NES), max(V(net)$NES), val_range/8), 
                             rev(brewer.pal(9,"PuRd")))
mycol = circlize::colorRamp2(seq(min(V(net)$NES), max(V(net)$NES), val_range/8), 
                             brewer.pal(9,"BuGn"))
V(net)$color= mycol(V(net)$NES)

#CHANGE THIS: if you have too many edges, clustering will not work + you'll end up with a hairball
net.sp = trim_edges(net, target_degree = .1)
clp <- cluster_leiden(net.sp, resolution_parameter = 2, 
                      objective_function = 'modularity',  n_iterations = 10)

#plots igraph network with clusters:
pdf(here(FIGDIR, 'plots','clustered.net.pdf'))
plot(clp, net.sp, vertex.label=NA)
dev.off()

# trim singletons and add the hub scores of each pathway
net_list = trim_clustered_nodes(net.sp, clp)
V(net_list$net)$hub_score = hub_score(net_list$net)$vector

## clean up the labels
vert.meta =  vertex_attr(net_list$net) %>% as.data.frame() 
to_label = vert.meta %>% clean_pathways()
to_label_num = vert.meta %>% dplyr::select(name, mem) %>% deframe()

## make 3 sets of plots
pdf(here(FIGDIR, 'plots','gsea_celltype_clustered_pathways.labeled.pdf'))
plotNetSimple(net_list$net, mark.group=net_list$vl, label.font = 2, label.cex = 2,
              vertex.label=to_label[V(net_list$net)$name], 
              vertex.shape = type_shape[V(net_list$net)$type])
dev.off()

pdf(here(FIGDIR, 'plots','gsea_celltype_clustered_pathways.numbered.pdf'))
plotNetSimple(net_list$net, mark.group=net_list$vl, 
              vertex.label=to_label_num[V(net_list$net)$name], 
              vertex.shape = type_shape[V(net_list$net)$type])
dev.off()

pdf(here(FIGDIR, 'plots','gsea_celltype_clustered_pathways.blank.pdf'))
plotNetSimple(net_list$net, mark.group=net_list$vl, 
              vertex.shape = type_shape[V(net_list$net)$type])
dev.off()

curenrich_clustered = curenrich %>% filter(OUD.v.CTL.in == 'All') %>%filter(NES > 0) %>% 
   mutate(cluster_number = to_label_num[pathway]) %>% arrange(cluster_number, padj) %>%
  relocate(cluster_number, pathway, description,.after = 'celltype') %>% 
  writexl::write_xlsx(here(FIGDIR, 'tables','gsea_celltype_clustered_pathways.numbered.xlsx'))


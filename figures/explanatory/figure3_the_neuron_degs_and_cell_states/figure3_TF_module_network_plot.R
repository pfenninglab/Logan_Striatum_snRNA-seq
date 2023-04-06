ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(rcartocolor)
library(data.table)
library(tidyverse)
library(here)
library(readxl)
library(wesanderson)
library(cividis) #  devtools::install_github("marcosci/cividis")
library(igraph)

DATADIR='data/tidy_data/differential_expression_analysis'

## make for this subdirs
PLOTDIR='figures/explanatory/figure3_the_neuron_degs_and_cell_states'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

###################################
# 0) pre-set plotting aesthetics
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 'All', 'Neuron', 'Glia',
               'Astrocytes', 'Endothelial', 'Microglia', 'Mural', 'Oligos', 
               'Oligos_Pre', 'Interneuron')
neuron_types = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome',
                 'D1/D2-Hybrid', 'Interneuron') %>% make.names()

in2mm<-25.4
my_theme = theme_classic(base_size = 5)
alpha = 0.05

sign_col = c('pink', 'lightgreen')
sign_col2 = c('darkred', 'darkgreen')

##############################################
# 1) read in the cell type OUD DEGs
rdasDir =file.path(DATADIR, 'rdas')
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds')
res = readRDS(save_res_fn)[c('Neuron', neuron_types)]

# number of cell types
length(res) 

# number of unique genes that are DEGs across any cell type
deg_df = res %>% rbindlist()
deg_signif_df = deg_df %>% filter(adj.P.Val.Between < alpha)
max_gene_expr = deg_df %>% arrange(desc(AveExpr)) %>% 
  distinct(gene, .keep_all = T) %>% dplyr::select(gene, AveExpr) %>% deframe()
deg_signif_sign = deg_signif_df %>% group_by(gene) %>% 
  summarise(sign = sign(mean(logFC))/2 + 1.5) %>% deframe()
deg_sign = deg_df %>% group_by(gene) %>% 
  summarise(sign = sign(mean(logFC))/2 + 1.5) %>% deframe()

#############################################
# 2) ead in the saved pseudo bulk TF regulon output 
stats_group  = here('figures/exploratory/AUCell_gene_set_activities', 'tables', 'OUD_Striatum_pyscenic_differential_TF_modules.byCelltype.xlsx') %>% 
  readxl::read_xlsx() %>% filter(celltype3 %in% neuron_types ) %>% 
  mutate(isSignif = p.adj < alpha )

signif_TF = stats_group %>% filter(isSignif ) %>% pull(TF) %>% unique()
signif_cell = stats_group %>% filter(isSignif ) %>% pull(celltype3) %>% unique()

## reorder the TFs by the difference in OUD
stats_group2 = stats_group %>% filter(TF %in% signif_TF, celltype3 %in% signif_cell) %>% 
  group_by(TF) %>% mutate(tmp = mean(statistic)) %>% ungroup() %>% 
  arrange(tmp) %>% mutate(TF = factor(TF, unique(TF)))

TF_sign= stats_group %>% filter(isSignif ) %>% group_by(TF) %>% 
  summarise(sign = sign(mean(estimate))/2 + 1.5) %>% deframe()
  
############################################
# 2) read in the TF-gene relationships
links = here('figures/exploratory/AUCell_gene_set_activities', 'tables', 'OUD_Striatum_pyscenic_detected_TF_modules.xlsx') %>% 
  readxl::read_xlsx() %>%
  dplyr::rename('TF' = 'V1', "gene" = "TargetGenes", 'edge.width' = 'importance') %>% 
  filter(TF %in% signif_TF) 

select_TF = c("NFATC2", 'PRRX1', 'BCLAF1', 'BCL11A', 'NR3C1', 'BACH1')

for(plot_TF in select_TF){
  set.seed(0)
  links2 = links %>% filter(TF %in% plot_TF) %>% 
    mutate(isDEG = gene %in% deg_signif_df$gene, 
           AveExpr = max_gene_expr[gene], 
           edge.width = edge.width**2) %>% 
    filter(edge.width > 1, AveExpr > 1) %>% 
    arrange(edge.width, AveExpr) %>% 
    top_n(50, edge.width)
  
  num_edges = nrow(links2)
  links2 %>% dplyr::count(isDEG)%>% print()
  
  nodes = links2 %>% dplyr::select(c(gene, TF, isDEG)) %>% 
    pivot_longer(cols = c(gene, TF), values_to = 'gene', names_to = 'type') %>% 
    mutate(type = type =='TF', 
           frame.color = ifelse(type, 'black', '#00000000'), 
           color = ifelse(type, sign_col2[TF_sign[gene]], 
                          ifelse(isDEG, sign_col2[deg_signif_sign[gene]], 
                                 sign_col[deg_sign[gene]])), 
           shape = 'rectangle', 
           size = 8*str_count(gene), 
           size2 = 15) %>% 
    relocate('gene', .before = everything()) %>% 
    arrange(desc(type)) %>% 
    distinct(gene, .keep_all = T)
  
  net <- graph_from_data_frame(d=links2, vertices=nodes, directed=T)
  V(net)$label.cex = .3
  V(net)$label.color = 'black'
  is_bipartite(net) #True
  
  ## make the plot of the gene modules
  plot_fn = here(PLOTDIR,'plots', paste0('figure3_TF-gene_modules_',plot_TF,'.pdf'))
  pdf(plot_fn, height = log10(num_edges)*30/in2mm, 
      width = log10(num_edges)*30/in2mm, onefile = T)
  par(mar=c(0,0,0,0))
  plot(net, edge.arrow.size=.1, edge.curved=.1, label.family = 'Helvetica', 
       layout=layout.fruchterman.reingold(net))
  dev.off()
}


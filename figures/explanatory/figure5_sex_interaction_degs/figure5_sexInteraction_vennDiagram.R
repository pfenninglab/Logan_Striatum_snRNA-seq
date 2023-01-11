ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggsci)
library(here)

library(ggvenn) # devtools::install_github("yanlinlin82/ggvenn")
library(VennDiagram)

DATADIR='data/tidy_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure5_sex_interaction_degs'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)
in2mm<-25.4

display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

###################################
# 0) pre-set plotting aesthetics
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes %>% make.names()

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'All', 'Neuron', 'Glia',
               'Astrocytes', 'Endothelial', 'Microglia', 'Mural', 'Oligos', 
               'Oligos_Pre', 'Interneuron')
othertypes_col = c(carto_pal(4, "Safe"), 
                   mypal = pal_npg('nrc')(3),
                   carto_pal(length(othertypes) -7 , "Vivid"))
names(othertypes_col) = othertypes
othertypes_col = othertypes_col[-c(1:4)]
typecolors = c(othertypes_col[1:3], subtypes_col, othertypes_col[c(10, 4:9)])

in2mm<-25.4
my_theme = theme_classic(base_size = 6)

#################################################
# 1) load in the DEG table from the big analyses
rdasDir =file.path(DATADIR,'differential_expression_analysis', 'rdas')

## DEGs of OUD vs. Control in either females or male
res_OUDwinSex = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds') %>% 
  readRDS() %>% lapply(function(df) df %>% mutate(group = 'OUDwinSex')) %>% rbindlist() %>% 
  pivot_longer(cols = -c(celltype, gene, group), names_sep = '_', names_to = c('metric', 'base_group'))

## DEGs of male vs. female in either OUD or control
res_SexWinOUD = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.SexWinOUD.rds') %>% 
  readRDS() %>% lapply(function(df) df %>% mutate(group = 'SexWinOUD')) %>% rbindlist()%>% 
  pivot_longer(cols = -c(celltype, gene, group), names_sep = '_', names_to = c('metric', 'base_group'))

## relabel the DEG comparisons to be meaningful
remap_groupings = c( OUD = 'M v F\nin OUD', CTL = 'M v F\nin CTL',
  SexF = 'OUD v CTL\nin F', SexM = 'OUD v CTL\nin M')
res = bind_rows(res_OUDwinSex, res_SexWinOUD) %>% 
  pivot_wider(names_from = 'metric', values_from = 'value') %>% 
  mutate(base_group = remap_groupings[base_group])



###############################################################
# 2) construct DEG list across 4 comparisons for venn diagram

## make a list for each groupping
alpha = 0.05
venn_list = split(res, f = res$celltype) %>% 
  lapply(function(x){
    x = x %>% filter(adj.P.Val.Between < alpha)
    split(x$gene, x$base_group)
  })


set_cols =  c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")
file_name = here(PLOTDIR, 'plots', paste0('figure5_sexInteraction_DEG_4-way_vennDiagram.pdf'))
pdf(file_name, height = 20/in2mm, width = 20/in2mm)
for (cell in names(venn_list)){
  file_name = here(PLOTDIR, 'plots', paste0('figure5_sexInteraction_DEG_4-way_vennDiagram.', cell, '.pdf'))
  pdf(file_name, height = 20/in2mm, width = 20/in2mm)
  display_venn( venn_list[[cell]],  fill = set_cols, main = cell, lty = 'blank',
                height = 20/in2mm, width = 20/in2mm, cex = .3, cat.cex = .25,  
                cat.fontface = "bold", resolution = 200, main.cex =.3)
  dev.off()
  display_venn( venn_list[[cell]],  fill = set_cols, main = cell, lty = 'blank',
                height = 20/in2mm, width = 20/in2mm, cex = .3, cat.cex = .25,  
                cat.fontface = "bold", resolution = 200, main.cex =.3)
}
dev.off()


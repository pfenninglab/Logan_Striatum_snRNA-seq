## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(rcartocolor)
library(ggupset)
library(svglite)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'
PLOTDIR='figures/exploratory/differential_expression_analysis/plots'

## load in the DEG table
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limmaRes_bigModel_N22.rds')
res = readRDS(save_res_fn)

alpha = 0.05
df = res %>% lapply(function(x){
  x[x$p_adj < alpha, ]
}) %>% rbindlist(idcol = 'celltype')%>% 
  filter(!celltype %in% c( "DSM.IV.OUD", "RegionPutamen", "SexM" ) )

###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural/Fibroblast', 'Oligos', 'Oligos_Pre')
othertypes = c('Interneuron',  'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = c(carto_pal(length(othertypes) , "Vivid"))
names(othertypes_col) = othertypes

typecolors = c(subtypes_col, othertypes_col)
names(typecolors) = make.names(names(typecolors))



###########################
# 1) make some upset plots

## neurons
plot_fn = here(PLOTDIR, 'OUD_Striatum_voomLimmaRes_N22.bycelltype3.upset.neuronDEG.pdf')
pdf(plot_fn, height = 4, width = 8)
df %>% filter(celltype %in% make.names(subtypes)) %>% 
  mutate(celltype = factor(celltype, names(typecolors))) %>% 
  group_by(gene) %>% summarize(celltype = list(celltype)) %>% 
  ggplot(aes(x = celltype)) + geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust='inward') +
  scale_x_upset() + 
  scale_y_continuous( name = "Overlapping DEGs") + 
  theme_bw(10) + theme(  plot.margin = margin(2,2,2,20, "mm"))
dev.off()


## glia
plot_fn = here(PLOTDIR, 'OUD_Striatum_voomLimmaRes_N22.bycelltype3.upset.gliaDEG.pdf')
pdf(plot_fn, height = 4, width = 8)
df %>% filter(celltype %in% make.names(othertypes)) %>% 
  filter(!grepl('Int', celltype) ) %>% 
  mutate(celltype = factor(celltype, names(typecolors))) %>% 
  group_by(gene) %>% summarize(celltype = list(celltype)) %>% 
  ggplot(aes(x = celltype)) + geom_bar() +
  geom_text(stat='count', aes(label=after_stat(count)), vjust='inward') +
  scale_x_upset() + 
  scale_y_continuous( name = "Overlapping DEGs") + 
  theme_bw(10) + theme(  plot.margin = margin(2,2,2,20, "mm"))
  dev.off()
  
  
  
  
  
  
  
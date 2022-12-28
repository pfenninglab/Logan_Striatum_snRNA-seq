ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(ggsci)
library(grid)
library(gridExtra) 

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## num parallel threads for Seurat
plan(sequential) # no parallel
options(future.globals.maxSize = 80e9)
library(here)

DATADIR='data/raw_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure2_the_celltype_annotations'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)


###################################
# 0) pre-set plotting aesthetics
in2mm<-25.4
my_theme = theme_classic(base_size = 6)

subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 'Astrocytes', 
               'Endothelial', 'Microglia', 'Mural', 'Oligos', 'Oligos_Pre', 
               'Interneuron')
othertypes_col = c(carto_pal(4, "Safe"), 
                   carto_pal(length(othertypes) -4 , "Vivid"))
names(othertypes_col) = othertypes
othertypes_col = othertypes_col[-c(1:4)]
othertypes_col = othertypes_col[c(7, 1:6)]
typecolors = c(subtypes_col, othertypes_col)

###################################
# 1) read in Logan snRNA dataset to plot
df = read_tsv(here('data/tidy_data/tables/BU_OUD_Striatum_refined_all_SeuratObj_N22.txt.gz'))
df = df %>% mutate(
  celltype3 = ifelse(grepl('Int', celltype3), 'Interneuron', celltype3),
  celltype3 = factor(celltype3, names(typecolors)))

df_summary = df %>% dplyr::select(ID:celltype3) %>% 
  group_by(ID, celltype3) %>%  mutate(numCell = n()) %>% ungroup()  %>% 
  distinct(ID, celltype3, .keep_all = T) %>% 
  group_by(ID) %>%  mutate(percent = numCell / sum(numCell) * 100) %>% ungroup()
  
df_summary %>% filter(ID =='C-1034') %>% dplyr::select(celltype3, percent)

df_prop = df_summary %>% group_by(celltype3) %>% 
  summarise(avg_percent = mean(percent)) %>% ungroup() %>% 
  mutate(avg_percent = avg_percent/sum(avg_percent)) %>% 
  write_tsv('data/tidy_data/tables/BU_OUD_Striatum_refined_celltype3_proportions.txt')


####################################################################
# 2) plot the diagonal matrix of opioid receptors and their endogenous peptide
sfig2_barplot = here(PLOTDIR, 'plots', 's2.3_celltype_proportion_barplots.pdf')
pdf(sfig2_barplot, width = 180/in2mm, height =  100/in2mm)
DotPlot( obj_merged, features = markerGenes, cols = c("lightgrey", 'limegreen'),
  cluster.idents = F, scale = T, scale.by = "radius") +
  my_theme + scale_y_discrete(limits = rev) +
  scale_x_discrete(position = "top") +
  guides(size = FALSE) +
  theme(legend.position = 'bottom', plot.title= element_blank(),
      axis.text.x = element_text(angle = -30, vjust = 0, hjust=1),
      axis.title = element_blank(), 
      legend.key.size = unit(2, "mm"),
      legend.spacing.x = unit(1, 'mm'))
dev.off()



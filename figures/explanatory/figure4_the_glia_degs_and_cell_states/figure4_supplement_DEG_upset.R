## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(data.table)
library(rcartocolor)
library(ComplexUpset) # devtools::install_github("krassowski/complex-upset")
library(cowplot)
library(ggsci)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


DATADIR='data/tidy_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure4_the_glia_degs_and_cell_states'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

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
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds')
res = readRDS(save_res_fn)

alpha = 0.05
df = res %>% lapply(function(x){
  x[x$adj.P.Val.Between < alpha, ]
}) %>% rbindlist()

###########################
# 2) make some upset plots

plot_fn = here(PLOTDIR, 'plots','s3.1_numDEG_bycelltype.upsetPlot.pdf')
pdf(plot_fn, height = 200/in2mm, width = 180/in2mm)

## neurons
ct1 = make.names(c('All', 'Neuron',subtypes,'Interneuron'))
dat1 = df %>% filter(celltype %in% ct1) %>% 
  dplyr::select(c('celltype', 'gene')) %>% mutate(val = 1) %>% 
  pivot_wider(names_from = celltype, values_from = val)%>% 
  mutate(across(-gene, ~replace_na(.x, 0))) %>%
  mutate(across(-gene, ~ifelse(.x==1, TRUE,FALSE)))
ct1 = ct1[ct1 %in% names(dat1)]

p1 = upset(dat1, rev(ct1), width_ratio=0.1, min_size=3, 
      wrap=TRUE, set_sizes=FALSE, sort_sets=FALSE) + 
  my_theme

## glia
ct2 = make.names(c('All', 'Glia', names(othertypes_col)[-c(1:3, 6, 10)]))
dat2 = df %>% filter(celltype %in% ct2) %>% 
  dplyr::select(c('celltype', 'gene')) %>% mutate(val = 1) %>% 
  pivot_wider(names_from = celltype, values_from = val)%>% 
  mutate(across(-gene, ~replace_na(.x, 0))) %>%
  mutate(across(-gene, ~ifelse(.x==1, TRUE,FALSE)))
ct2 = ct2[ct2 %in% names(dat2)]

p2 = upset(dat2, rev(ct2), width_ratio=0.1, min_size=3, 
           wrap=TRUE, set_sizes=FALSE, sort_sets=FALSE) + my_theme

## put it together
p1 /   p2
dev.off()


#####################################################################
# 3) save the list of DEGs that are shared across multiple cell types
df_neuron = df %>% filter(celltype %in% ct1)
df_glia = df %>% filter(celltype %in% ct2)
df_all = df %>% arrange(adj.P.Val.Between)

df_list = list('All_overlap' = df_all,
               'Neuron_overlap' = df_neuron, 
               'Glia_overlap' = df_glia) %>% 
  lapply(function(df){ 
   df %>% 
      arrange(adj.P.Val.Between) %>% group_by(gene) %>% 
      mutate(
        numCelltype = n(),
        isSameFCdir = var(sign(logFC))==0,
        logFC = mean(logFC),
        celltype = paste(celltype, collapse = ', ')
      ) %>% 
      dplyr::relocate(celltype, .after = isSameFCdir) %>% 
      top_n(1, wt = -adj.P.Val.Between) %>% ungroup() %>% 
      filter(numCelltype > 1) %>% 
      arrange(desc(numCelltype), desc(abs(logFC)))
  })

save_overlap = here(PLOTDIR, 'tables','s3.2_table_numDEG_overlap_bycelltype.xlsx')
df_list %>% writexl::write_xlsx(save_overlap)


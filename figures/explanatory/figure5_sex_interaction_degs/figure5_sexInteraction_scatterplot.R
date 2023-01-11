ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggsci)
library(here)

DATADIR='data/tidy_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure5_sex_interaction_degs'
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
res_cell = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% 
  readRDS() %>% lapply(function(df) df %>% mutate(group = 'celltype')) %>% rbindlist()
res_sex = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.sexInteraction.rds') %>% 
  readRDS() %>% lapply(function(df) df %>% mutate(group = 'sex')) %>% rbindlist()

res = bind_rows(res_cell, res_sex)

######################################
## 2) count DEGs within each grouping
alpha = 0.05
df = res %>% group_by(gene, celltype) %>% 
  filter(any(adj.P.Val.Between< alpha)) %>% 
  ungroup() 

df_wide= df %>% dplyr::select(celltype:t, adj.P.Val.Between, group) %>% 
  pivot_wider(names_from='group', values_from = c(logFC,t, AveExpr, adj.P.Val.Between)) %>% 
  mutate(
    celltype = factor(celltype, names(typecolors)),
    signif_group = case_when(
    adj.P.Val.Between_celltype < alpha & adj.P.Val.Between_sex<alpha ~ 'signif.both', 
    adj.P.Val.Between_celltype < alpha ~ 'signif.celltype', 
    adj.P.Val.Between_sex < alpha ~ 'signif.sex'))

df_wide = df_wide %>% 
  mutate(logFC_diff = logFC_sex - logFC_celltype, 
         logAdjP_diff = -log10(adj.P.Val.Between_sex) + log10(adj.P.Val.Between_celltype),
         logAdjP_diff = pmax(logAdjP_diff, 0) ) %>% 
  arrange(desc(logAdjP_diff), desc(logFC_diff))


df_wide %>% filter(celltype %in% c('Microglia', 'Oligo')) %>% as.data.frame() %>% head()
##########################################
## 2) make the scatter plot of the log2FC's
plot_fn = here(PLOTDIR,'plots', 'sX.2_mumDEG_sexInteraction_scatterplot.logFC.pdf')
pdf(plot_fn, height = 220/in2mm, width = 180/in2mm)
ggplot(df_wide, aes(x = logFC_celltype, y = logFC_sex))+
  geom_point(pch = 21, aes(fill = signif_group)) + 
  geom_hline(yintercept = 0, color = 'red') +
  geom_vline(xintercept = 0, color = 'red') +
  facet_wrap( ~ celltype, scales = 'fixed')+
  xlab('Log2FC by OUD Dx in each cell type, (OUD-CTL)') + 
  ylab('Log2FC interaction between Sex and OUD Dx, (OUD.M-OUD.F) - (CTL.M - CTL.F)')+
  scale_y_continuous(labels = abs) + # so negative sign doesn't show
  my_theme + theme(legend.position = 'bottom')
dev.off()


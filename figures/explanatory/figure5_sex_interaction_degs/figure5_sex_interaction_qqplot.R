ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(tidymodels)
library(broom)
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

in2mm<-25.4; alpha = 0.05
my_theme = theme_classic(base_size = 6)

#################################################
# 1) load in the DEG table from the big analyses
rdasDir =file.path(DATADIR,'differential_expression_analysis', 'rdas')

## DEGs of OUD vs. Control in either females or male
res = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.SexSpecific.rds') %>% 
  readRDS() %>% lapply(function(df) df %>% mutate(group = 'OUDwinSex')) %>% 
  rbindlist() 

## calculate the average number of interaction effects
res %>% filter(adj.P.Val.Between_Interaction < alpha) %>% 
  filter(!celltype %in% c('Neuron', 'All', 'Glia')) %>% 
  mutate(cell_class = ifelse(celltype %in% c(make.names(subtypes), 'Interneuron'), 'Neuron', 'Glia')) %>%
  group_by(cell_class, celltype) %>% summarise(n = n()) %>%
  group_by(cell_class) %>% 
  summarise(mean = mean(n), sem = sd(n)/sqrt(n()))


##########################################
## 2) make the scatter plot of the log2FC's
typecolors = typecolors[names(typecolors) %in% res2$celltype]

## calculate the observed and expected p-Values for qqplot
res2 = res %>% 
  arrange(adj.P.Val.Between_Interaction) %>% 
  group_by(celltype) %>%
  mutate(  observed = -log10(adj.P.Val.Between_Interaction),
           expected = -log10(ppoints(length(adj.P.Val.Between_Interaction)))) %>% 
  ungroup()

plot_fn = here(PLOTDIR,'plots', 'sX.2_mumDEG_sexInteraction.qqplot.pdf')
pdf(plot_fn, height = 121/in2mm, width = 121/in2mm)
ggplot(res2, aes(x = expected, y = observed))+
  geom_abline(intercept = 0, slope = 1, alpha = 0.5, lty = 'dashed')+
  geom_line(aes(color = celltype)) + 
  # geom_point(pch = 21, aes(fill = celltype)) + 
  scale_fill_manual(values = typecolors, guide = 'none') +
  scale_color_manual(values = typecolors, guide = 'none') +
  facet_wrap( ~ celltype, scales = 'fixed')+
  ylab('Observed -log10(P-Values)') +
  xlab('Expected -log10(P-Values)') + 
  my_theme + theme(legend.position = 'right', 
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 10, hjust=1))
dev.off()







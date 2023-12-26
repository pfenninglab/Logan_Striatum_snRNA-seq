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
  rbindlist() %>% filter(!celltype %in% c('All', 'Neuron', 'Glia'))

df_overlap = res %>% 
  mutate(
    cell_class = case_when(
      celltype %in% c(make.names(subtypes), 'Interneuron') ~ 'Neuron', 
      celltype %in% make.names(othertypes) ~ 'Glia'
    ),
    signif_group = case_when(
    adj.P.Val.Between_SexF < alpha & adj.P.Val.Between_SexM < alpha ~ 'in Both', 
    adj.P.Val.Between_SexF < alpha ~ 'in Females only', 
    adj.P.Val.Between_SexM < alpha ~ 'in Males only'
  )) %>% filter(!is.na(signif_group)) %>% 
  group_by(cell_class, celltype) %>% dplyr::count(signif_group) %>% 
  mutate(
    cell_class = factor(cell_class, levels = c('Neuron', 'Glia')),
    celltype = factor(celltype, levels = make.names(c(subtypes, othertypes))),
  )
df_overlap %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'sX.2_mumDEG_sexInteraction_venn.sourceData.xlsx'))

df_overlap %>% group_by(cell_class, signif_group) %>% 
  summarise(mean = mean(n), sem = sd(n)/sqrt(n()))

##########################################
## 2) make the scatter plot of the log2FC's
typecolors = typecolors[names(typecolors) %in% df_overlap$celltype]

plot_fn = here(PLOTDIR,'plots', 'sX.2_mumDEG_sexInteraction_venn.barplot.pdf')
pdf(plot_fn, height = 50/in2mm, width = 121/in2mm)
ggplot(df_overlap, aes(x = celltype, y = n, fill = celltype))+
  geom_bar(stat = 'identity', position = position_dodge()) + 
  geom_text(aes(label = n, color = signif_group), color = 'black',
            vjust = 'inward', position = position_dodge(.9), size = 2) + 
  scale_fill_manual(values = typecolors, guide = 'none') +
  facet_grid(signif_group ~ cell_class, scales = 'free_x', space = 'free')+
  ylab('# of OUD vs. UC DEGs by Sex (FDR < 0.05)')+
  scale_y_continuous(labels = abs) + # so negative sign doesn't show
  my_theme + theme(legend.position = 'right', 
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 10, hjust=1))
dev.off()



#######################################
## 3) count the interaction effect DEGs 

df_overlap2 = res %>% 
  filter(adj.P.Val.Between_Interaction < alpha) %>% 
  mutate(
    cell_class = case_when(
      celltype %in% c(make.names(subtypes), 'Interneuron') ~ 'Neuron', 
      celltype %in% make.names(othertypes) ~ 'Glia'
    ),
    bias = case_when(
      abs(dir_SexM) > abs(dir_SexF) ~ 'M',
      abs(dir_SexM) < abs(dir_SexF) ~ 'F',
    ),
    signif_group = case_when(
      logFC_Interaction > 0 & bias =='M' ~ 'Up in OUD\nin Males', 
      logFC_Interaction < 0 & bias =='F' ~ 'Down in OUD\nin Females', 
      logFC_Interaction > 0 & bias =='F' ~ 'Up in OUD\nin Females', 
      logFC_Interaction < 0 & bias =='M' ~ 'Down in OUD\nin Males', 
    ), 
    signif_group = factor(signif_group, levels = 
                            c('Up in OUD\nin Females', 'Down in OUD\nin Females', 'Up in OUD\nin Males', 'Down in OUD\nin Males'))
    ) %>% 
  filter(!is.na(signif_group)) %>% 
  group_by(cell_class, celltype) %>% count(signif_group) %>% 
  mutate(
    cell_class = factor(cell_class, levels = c('Neuron', 'Glia')),
    celltype = factor(celltype, levels = make.names(c(subtypes, othertypes))),
  )

plot_fn = here(PLOTDIR,'plots', 'sX.8_mumDEG_sexInteraction.barplot.pdf')
pdf(plot_fn, height = 80/in2mm, width = 121/in2mm)
ggplot(df_overlap2, aes(x = celltype, y = n, fill = celltype))+
  geom_bar(stat = 'identity', position = position_dodge()) + 
  geom_text(aes(label = n, color = signif_group), color = 'black',
            vjust = 'inward', position = position_dodge(.9), size = 2) + 
  scale_fill_manual(values = typecolors, guide = 'none') +
  facet_grid(signif_group ~ cell_class, scales = 'free_x', space = 'free')+
  ylab('# of OUD vs. UC DEGs with Sex Interaction (FDR < 0.05)')+
  scale_y_continuous(labels = abs) + # so negative sign doesn't show
  my_theme + theme(legend.position = 'right', 
                   axis.title.x = element_blank(),
                   axis.text.x = element_text(angle = 10, hjust=1))
dev.off()


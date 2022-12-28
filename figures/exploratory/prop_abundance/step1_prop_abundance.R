## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

library(tidymodels)
library(broom.mixed)
library(ggsci)
library(lme4)
library(lmerTest)
library(qvalue)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


## make for this subdirs
PLOTDIR='figures/exploratory/prop_abundance'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

###################################
# 0) pre-set plotting aesthetics
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

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

dx_col = setNames(pal_npg(palette = c("nrc"), alpha = 1)(2), c('CTL', 'OUD'))

in2mm<-25.4
my_theme = theme_classic(base_size = 6)

##################################################
# 1) load in full dataset cell type labels for plot

## load the unfiltered QC table
meta_df = here('data/tidy_data/tables',
               paste0("BU_OUD_Striatum_refined_all_SeuratObj_N22.txt.gz")) %>%
  read_tsv(show_col_types = FALSE) 

pd = meta_df %>% dplyr::select(ID:DSM.IV.CUD) %>% distinct(ID, .keep_all= TRUE)

prop_df = meta_df %>% 
  mutate(celltype3 = ifelse(grepl('Int', celltype3), 'Interneuron', celltype3),
         celltype3 = factor(celltype3, names(typecolors)) %>% droplevels()) %>% 
  group_by(ID, celltype3) %>% 
  summarise(numCelltype = n()) %>% group_by(ID) %>% 
  mutate(propCelltype = numCelltype / sum(numCelltype) * 100) %>% 
  distinct(.keep_all = TRUE) %>% 
  ungroup() %>% left_join(pd) 

avg_df = prop_df %>% group_by(celltype3) %>% summarise(avg_prop = mean(propCelltype))

## test if there's a difference in cell type proportions w/ OUD
lm_prop = prop_df %>% 
  nest(data = -c(celltype3)) %>% 
  mutate(
    test = map(data, ~ lmer(propCelltype ~ DSM.IV.OUD + Region +
                              Age + Sex + PMI + RIN + (1|Case), data = .x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% select(-data, -test) %>% 
  filter(grepl('OUD', term)) %>% 
  mutate(fdr = p.adjust(p.value, 'fdr')) %>% 
  inner_join(avg_df) %>% arrange(p.value)

lm_prop %>% as.data.frame() %>% head()


out_fn = here(PLOTDIR, 'tables', 'cell_type_proportion_byOUD.boxplot.xlsx')
lm_prop %>% writexl::write_xlsx(out_fn)


## make the plots
pdf(here(PLOTDIR, 'plots', 'cell_type_proportion_byOUD.boxplot.pdf'), width = 4.5, height = 2)
ggplot(prop_df, aes(x = DSM.IV.OUD, y = propCelltype)) +
  geom_boxplot(outlier.shape = NA, aes(color = DSM.IV.OUD)) + 
  geom_point(alpha = 0.8, aes( fill = DSM.IV.OUD, shape = Sex),
             position = position_jitterdodge(jitter.width = .25), size = 1) + 
  geom_text(data = lm_prop, aes(label = label, x = 1.5, y = avg_prop*1.5), 
            hjust="inward") +
  scale_fill_manual(values =c('CTL' = 'gray', 'OUD' = 'black')) +
  scale_shape_manual(values =c('F' = 21, 'M' = 22)) +
  scale_color_manual(values =c('CTL' = 'gray', 'OUD' = 'black')) +
  facet_wrap(~celltype3, scale = 'free_y', nrow = 2)+
  theme_bw(base_size = 5) +
  ylab('Percent of Sample') + 
  scale_y_continuous(expand = expansion(mult = c(.1, .3)))+
  theme(legend.position = 'bottom',
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()




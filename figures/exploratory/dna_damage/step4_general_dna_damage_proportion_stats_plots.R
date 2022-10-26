## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

## statistical analyses
library(tidymodels)
library(lme4)
library(broom.mixed)
library(lmerTest)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

PLOTDIR='figures/exploratory/dna_damage'

#########################################
# 1) load in the DNA damage estimates and annotations
save_file = here(PLOTDIR, 'rdas', 
                 'OUD_Striatum_refined_all_SeuratObj_N22.meta.DNAdam2.rds')
celltype_df = readRDS(save_file)

########################
# 2) average per sample
dam_per_sample = celltype_df%>% 
  mutate(Case = ss(ID, '-', 2)) %>% group_by(Case) %>% 
  mutate(numCell = n(), DNA_dam_prop_all = sum(DNA_dam_class_all =='damaged') / numCell) %>% 
  ungroup() %>% distinct(Case, .keep_all = T) 

## stats, no difference in proportion of damaged cells
modBySample = glm(DNA_dam_prop_all ~ DSM.IV.OUD+ Age + Sex + PMI + RIN + numCell,
                 data = dam_per_sample, family = 'binomial')
summary(modBySample) 

##############################################
# 3) stratified by cell type and sample
dam_per_cellxSample = celltype_df %>% 
  arrange(celltype3) %>% 
  mutate(
    celltype4 = as.character(celltype3),
    celltype4 = ifelse(grepl('D1-M|D1-S', celltype4), 'D1-MSN', celltype4),
    celltype4 = ifelse(grepl('D2-M|D2-S', celltype4), 'D2-MSN', celltype4),
    celltype4 = ifelse(grepl('Int-', celltype4), 'Interneurons', celltype4)) %>% 
  group_by(Case, Region, celltype4) %>% 
  mutate_if(is.numeric, mean) %>% 
  mutate(numCell = n(), 
         damaged = sum(DNA_dam_class_all =='damaged'), 
         undamaged = sum(DNA_dam_class_all =='undamaged'), 
         DNA_dam_prop_all = damaged / numCell) %>% 
  ungroup() %>% distinct(Case, Region, celltype4, .keep_all = T) %>% 
  group_by(celltype4) %>% mutate(tmp = mean(DNA_dam_val_all)) %>% 
  ungroup() %>% arrange(desc(DNA_dam_val_all)) %>%   
  mutate(celltype4 = factor(celltype4, unique(celltype4)))

## for plotting later, get average proportion of cell types
dam_per_cellxSample2 = dam_per_cellxSample %>% 
  mutate(damaged= damaged / numCell * 100, 
         undamaged= undamaged / numCell * 100) %>% 
  group_by(DSM.IV.OUD, celltype4) %>% 
  summarise_if(is.numeric, mean) %>% 
  ungroup() %>% pivot_longer(cols = c('damaged', 'undamaged')) %>% 
  mutate(name = factor(name, c('undamaged', 'damaged')),
         annotation = paste0(sprintf("%2.1f", value), "%"))

## filter out cell types without any damaged cells by AUCell cutoff
dam_per_cellxSample = dam_per_cellxSample  %>% 
  ## z-normalize to rescale these numeric values for regression
  mutate_at(all_of(c('Age', 'PMI', 'RIN', 'numCell')), 
            ~ (. - mean(.))/sd(.)) %>% 
  # arrange(!grepl('D1', celltype4), !grepl('D2', celltype4), 
  #         !grepl('Int', celltype4), celltype4) %>% 
  group_by(celltype4, Region) %>% 
  mutate(toKeep = mean(DNA_dam_prop_all[DSM.IV.OUD=='OUD']) != 
                         mean(DNA_dam_prop_all[DSM.IV.OUD=='CTL'])) %>% 
  ungroup() %>% filter(toKeep)

## stats using mixed effect averages of cell type DNA dam scores per celltype, region, and patient
lmer(DNA_dam_prop_all ~ celltype4 + celltype4:DSM.IV.OUD+ Age + Sex + 
       PMI + RIN + numCell + Region + (1|Case), 
     data = dam_per_cellxSample) %>% summary()

modBySampleAndCell = lmer(DNA_dam_prop_all ~ celltype4 + celltype4:DSM.IV.OUD+ Age + Sex + 
                            PMI + RIN + numCell + Region + (1|Case), 
                          data = dam_per_cellxSample) %>% tidy() %>% 
  filter(grepl('^celltype4', term) & grepl('OUD$', term)) %>% arrange(p.value) %>% 
  dplyr::select(-c('effect', 'group')) %>% 
  mutate(celltype = ss(term, '4|:', 2))
modBySampleAndCell %>% as.data.frame() %>% relocate('p.value', .after = 'term')


################################################################
# 4) plot all the cells together general DNA damage proportions

## plots all the subtypes
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamAllProp.perSampleByCell.pdf'), width = 4.5, height = 1)
ggplot(dam_per_cellxSample2 %>% filter(name == 'damaged'), 
       aes(x =DSM.IV.OUD, y = value, fill = DSM.IV.OUD)) +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = annotation, y = value + 1), size = 1.5)+
  scale_fill_manual(values =c('gray', '#2162B0'), guide = 'none') +
  facet_wrap(~celltype4, scale = 'fixed', nrow = 1)+ 
  # coord_cartesian(ylim = c(0, .15))+
  theme_bw(base_size = 5) +
  ylab('General Damaged %') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()



###################################################
# 5) split the proportion plots by glia and neurons

## get split
glia = dam_per_cellxSample2 %>% filter(!duplicated(celltype4)) %>% 
  filter(grepl('As|lia|lig|Mur', celltype4)) %>% pull(celltype4)
neur = dam_per_cellxSample2 %>% filter(!duplicated(celltype4)) %>% 
  filter(!grepl('As|lia|lig|Mur', celltype4)) %>% pull(celltype4)

## plot the glia
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamAllProp.perSampleByGlia.pdf'), width = 2.5, height = 1)
ggplot(dam_per_cellxSample2 %>% filter(name == 'damaged' & celltype4 %in% glia), 
       aes(x =DSM.IV.OUD, y = value, fill = DSM.IV.OUD)) +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = annotation, y = value + 1), size = 1.5)+
  scale_fill_manual(values =c('gray', '#2162B0'), guide = 'none') +
  facet_wrap(~celltype4, scale = 'fixed', nrow = 1)+ 
  # coord_cartesian(ylim = c(0, .15))+
  theme_bw(base_size = 5) +
  ylab('General Damaged %') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()


## plot the neurons
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamAllProp.perSampleByNeuron.pdf'), width = 2, height = 1)
ggplot(dam_per_cellxSample2 %>% filter(name == 'damaged' & celltype4 %in% neur), 
       aes(x =DSM.IV.OUD, y = value, fill = DSM.IV.OUD)) +
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = annotation, y = value + 1), size = 1.5)+
  scale_fill_manual(values =c('gray', '#2162B0'), guide = 'none') +
  facet_wrap(~celltype4, scale = 'fixed', nrow = 1)+ 
  # coord_cartesian(ylim = c(0, .15))+
  theme_bw(base_size = 5) +
  ylab('General Damaged %') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()


###########################
# 5) export stats to table 
list('DNA_dam_prop_allortion_by_sample' = modBySample %>% tidy(), 
     'DNA_dam_prop_allortion_by_sampleCelltype' = modBySampleAndCell) %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_dnaDamAllProp.perSampleByCell.xlsx'))


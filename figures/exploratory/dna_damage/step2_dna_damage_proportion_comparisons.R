## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## AUCell object for DNA damage estimates
library(AUCell)

## for plots
library(Nebulosa)
library(ggbeeswarm)

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
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F)
dir.create(here(PLOTDIR, 'rdas'), showWarnings = F)

###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH', 'Int-SST', 'Int-TH', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = carto_pal(length(othertypes), "Vivid")
names(othertypes_col) = othertypes


#########################################
# 1) load in the DNA damage estimates and annotations
save_file = here(PLOTDIR, 'rdas', 
                 'BU_OUD_Striatum_refined_all_SeuratObj_N22.meta.DNAdam.rds')
celltype_df = readRDS(save_file)

########################
# 2) average per sample
dam_per_sample = celltype_df%>% 
  mutate(Case = ss(ID, '-', 2)) %>% group_by(Case) %>% 
  mutate(numCell = n(), DNA_dam_prop = sum(DNA_dam_class =='damaged') / numCell) %>% 
  ungroup() %>% distinct(Case, .keep_all = T) 

## stats, no difference in proportion of damaged cells
modBySample = glm(DNA_dam_prop ~ DSM.IV.OUD+ Age + Sex + PMI + RIN + numCell,
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
         damaged = sum(DNA_dam_class =='damaged'), 
         undamaged = sum(DNA_dam_class =='undamaged'), 
         DNA_dam_prop = damaged / numCell) %>% 
  ungroup() %>% distinct(Case, Region, celltype4, .keep_all = T) %>% 
  ## z-normalize to rescale these numeric values for regression
  mutate_at(all_of(c('Age', 'PMI', 'RIN', 'numCell')), 
            ~ (. - mean(.))/sd(.)) %>% 
  # arrange(!grepl('D1', celltype4), !grepl('D2', celltype4), 
  #         !grepl('Int', celltype4), celltype4) %>% 
  group_by(celltype4) %>% mutate(tmp = mean(DNA_dam_val)) %>% 
  ungroup() %>% arrange(desc(DNA_dam_val)) %>%   
  mutate(celltype4 = factor(celltype4, unique(celltype4)))

dam_per_cellxSample2=dam_per_cellxSample %>%  pivot_longer(cols = c('damaged', 'undamaged'))

## filter out cell types without any damaged cells by AUCell cutoff
dam_per_cellxSample = dam_per_cellxSample %>% 
  group_by(celltype4, Region) %>% 
  mutate(toKeep = mean(DNA_dam_prop[DSM.IV.OUD=='OUD']) != 
                         mean(DNA_dam_prop[DSM.IV.OUD=='CTL'])) %>% 
  ungroup() %>% filter(toKeep)

## stats using mixed effect averages of cell type DNA dam scores per celltype, region, and patient
lmer(DNA_dam_prop ~ celltype4 + celltype4:DSM.IV.OUD+ Age + Sex + 
       PMI + RIN + numCell + Region + (1|Case), 
     data = dam_per_cellxSample) %>% summary()

modBySampleAndCell = lmer(DNA_dam_prop ~ celltype4 + celltype4:DSM.IV.OUD+ Age + Sex + 
                            PMI + RIN + numCell + Region + (1|Case), 
                          data = dam_per_cellxSample) %>% tidy() %>% 
  filter(grepl('^celltype4', term) & grepl('OUD$', term)) %>% arrange(p.value) %>% 
  dplyr::select(-c('effect', 'group')) %>% 
  mutate(celltype = ss(term, '4|:', 2))
modBySampleAndCell %>% as.data.frame() %>% relocate('p.value', .after = 'term')

## plots
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_dnaDamProp.perSampleByCell.pdf'), width = 4.5, height = 1)
ggplot(dam_per_cellxSample2, 
       aes(x =DSM.IV.OUD, y = value, fill = name)) +
  geom_bar(position="fill", stat = 'identity') + 
  # geom_point(pch = 21, alpha = 0.5, position = position_jitter(width = 0.3)) + 
  scale_fill_manual(values =c('undamaged' = 'gray', 'damaged' = 'black')) +
  facet_wrap(~celltype4, scale = 'fixed', nrow = 1)+ 
  coord_cartesian(ylim = c(0.85, 1))+
  theme_bw(base_size = 5) +
  ylab('DNA Damaged Proportion') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()


###########################
# 5) export stats to table 
list('DNA_dam_proportion_by_sample' = modBySample %>% tidy(), 
     'DNA_dam_proportion_by_sampleCelltype' = modBySampleAndCell) %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'BU_OUD_Striatum_dnaDamProp.perSampleByCell.xlsx'))


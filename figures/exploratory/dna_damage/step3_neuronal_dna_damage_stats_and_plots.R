## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## for plot
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
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F)
dir.create(here(PLOTDIR, 'rdas'), showWarnings = F)

###############################################################
# 1) load in the DNA damage estimates for glia from the general 
# DNA dam vs. undamaged

save_file = here(PLOTDIR, 'rdas', 'BU_OUD_Striatum_refined_all_SeuratObj_N22.meta.DNAdam2.rds')
celltype_df = readRDS(save_file) %>% filter(!is.na(DNA_dam_val_neu)) %>% 
  filter(grepl('^D[1-2]|Int', celltype3))

celltype_df %>% dplyr::count(celltype3)

########################
# 2) average per sample
dam_per_sample = celltype_df %>% 
  mutate(Case = ss(ID, '-', 2)) %>% group_by(Case) %>% 
  mutate(numCell = n(), DNA_dam_prop = sum(DNA_dam_class_neu =='damaged') / numCell) %>% 
  mutate_if(is.numeric, mean) %>% 
  ungroup() %>% distinct(Case, .keep_all = T) 

dam_per_sample %>% dplyr::select(Case, DNA_dam_val_neu, DSM.IV.OUD, 
                                 Age, Sex, PMI, RIN, numCell) %>% 
  writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_dnaDamNeuVal.perSample.soureData.xlsx'))

## stats no difference in cortical neuronal DNA damage signature in MSNs
modBySample = lm(DNA_dam_val_neu ~ DSM.IV.OUD+ Age + Sex + PMI + RIN + numCell,
                 data = dam_per_sample)
summary(modBySample) 

## plots
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamNeuVal.perSample.pdf'), width = .75, height = 1)
ggplot(dam_per_sample, aes(x =DSM.IV.OUD, y = DNA_dam_val_neu)) +
  geom_boxplot(outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
  geom_line(alpha = 0.9, aes(x =DSM.IV.OUD, y = DNA_dam_val_neu, group=Pair)) + 
  geom_point(pch = 20, alpha = 0.9, aes(x =DSM.IV.OUD, y = DNA_dam_val_neu, color = Sex)) + 
  scale_color_manual(values =c('chocolate', 'aquamarine')) +
  scale_fill_manual(values =c('white', '#C9178D95'), guide = 'none') +
  theme_bw(base_size = 5) +
  ylab('Neuronal Damage Score') + xlab('')+
  theme(legend.position = 'top', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()


##############################################
# 3) compute group average DNA damage ratios

## stratified by cell type and sample
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
         DNA_dam_prop = sum(DNA_dam_class_neu =='damaged') / numCell) %>% 
  ungroup() %>% distinct(Case, Region, celltype4, .keep_all = T) %>% 
  ## z-normalize to rescale these numeric values for regression
  mutate_at(all_of(c('Age', 'PMI', 'RIN', 'numCell')), 
            ~ (. - mean(.))/sd(.)) %>% 
  # arrange(!grepl('D1', celltype4), !grepl('D2', celltype4), 
  #         !grepl('Int', celltype4), celltype4) %>% 
  group_by(celltype4) %>% mutate(DNA_dam_val_neu_avg = mean(DNA_dam_val_neu)) %>% 
  ungroup() %>% arrange(desc(DNA_dam_val_neu_avg)) %>%   
  mutate(celltype4 = factor(celltype4, unique(celltype4)))

dam_per_cellxSample %>% dplyr::select(Case, celltype4, DNA_dam_val_neu, DSM.IV.OUD, 
                                 Age, Sex, PMI, RIN, numCell) %>% 
  writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_dnaDamNeuVal.perSampleByCell.soureData.xlsx'))


## take a look
dam_per_cellxSample %>% count(Region, Case)
dam_per_cellxSample$DNA_dam_prop %>% summary()

## stats using mixed effect averages of cell type DNA dam scores per celltype, region, and patient
lmer(DNA_dam_val_neu ~ celltype4 + celltype4:DSM.IV.OUD+ Age + Sex + 
       PMI + RIN + numCell + Region + (1|Case), 
     data = dam_per_cellxSample) %>% summary()

# Pr(>|t|)    
# (Intercept)                         1.83e-08 ***
#   celltype4Interneurons:DSM.IV.OUDOUD  0.03208 *  
#   celltype4D2-MSN:DSM.IV.OUDOUD        0.49600    
#   celltype4D1/D2-Hybrid:DSM.IV.OUDOUD  0.49528    
#   celltype4D1-MSN:DSM.IV.OUDOUD        0.25020

modBySampleAndCell = lmer(DNA_dam_val_neu ~ celltype4 + celltype4:DSM.IV.OUD+ Age + Sex + 
                            PMI + RIN + numCell + Region + (1|Case), 
                          data = dam_per_cellxSample) %>% tidy() %>% 
  filter(grepl('^celltype4', term) & grepl('OUD$', term)) %>% arrange(p.value) %>% 
  dplyr::select(-c('effect', 'group')) %>% 
  mutate(celltype4 = ss(term, '4|:', 2) %>% factor(levels(dam_per_cellxSample$celltype4)), 
         annotation = case_when(
           p.value < 0.001 ~ '***',
           p.value < 0.01 ~ '**',
           p.value < 0.05 ~ '*',
           TRUE ~ "")) %>% 
  inner_join(dam_per_cellxSample %>% filter(!duplicated(celltype4)) %>% 
               rename('y_position' = 'DNA_dam_val_neu_avg') %>% 
               dplyr::select(celltype4, y_position))

modBySampleAndCell %>% as.data.frame() %>% relocate('p.value', .after = 'term') 

##############################################
# 4) compute group average DNA damage ratios

## plots all the cells
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamNeuVal.perSampleByCell.pdf'), width = 2, height = 1)
ggplot(dam_per_cellxSample, 
       aes(x =DSM.IV.OUD, y = DNA_dam_val_neu)) +
  geom_violin(size = .25, aes(fill = DSM.IV.OUD)) + 
  geom_boxplot(width = 0.4, size = .25, outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
  scale_fill_manual(values =c('white', '#C9178D95'), guide = 'none') +
  geom_text(data = modBySampleAndCell, size = 2, x = 1.5,
              aes(y = y_position + .01, label = annotation)) +
  facet_wrap(~celltype4, scale = 'fixed', nrow = 1)+
  theme_bw(base_size = 5) +
  ylab('Neuronal Damage Score') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()

###########################
# 5) export stats to table 
list('DNA_dam_score_by_sample' = summary(modBySample) %>% tidy(), 
     'DNA_dam_score_by_sampleCelltype' = modBySampleAndCell) %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_dnaDamNeuVal.perSampleByCell.xlsx'))


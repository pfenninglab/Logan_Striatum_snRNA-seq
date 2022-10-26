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
# 1) load in the DNA damage estimates
dnaDam_val = here('data/tidy_data/Seurat_projects/AUCell', 
                  '22_07_07_cells_AUC.rds') %>% readRDS() %>% 
  getAUC() %>% as.matrix() %>% t() %>% as.data.frame() 

dnaDam_class = here('data/tidy_data/Seurat_projects/AUCell', 
                    '22_07_07_cells_assignment.rds') %>% 
  readRDS() %>% unlist()


#####################################
# 2) load in cell type annotations 
obj = here('data/tidy_data/Seurat_projects', 
           "OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

## ordered cell types
celltypes = obj[[]] %>% filter(!duplicated(celltype3)) %>% 
  mutate(class = case_when(grepl('^D[1-2]', celltype3) ~ 'MSN', 
                           grepl('^Int', celltype3) ~ 'INT', 
                           TRUE ~ 'GLIA'), 
         class = factor(class, c('MSN', 'INT', 'GLIA'))) %>% 
  arrange(class, celltype3) %>% pull(celltype3)

## load the unfiltered QC table
celltype_df = obj[[]] %>% dplyr::select(-contains('integrated_snn_res')) %>% 
  mutate(celltype3 = factor(celltype3, celltypes))

celltype_df = cbind(celltype_df, 'DNA_dam_val' = dnaDam_val[rownames(celltype_df),]) %>% 
  mutate(DNA_dam_class = ifelse(rownames(celltype_df) %in% dnaDam_class, 'damaged', 'undamaged')) %>% filter(!is.na(DNA_dam_val))
  
save_file = here(PLOTDIR, 'rdas', 'OUD_Striatum_refined_all_SeuratObj_N22.meta.DNAdam.rds')
celltype_df %>% saveRDS(save_file)

cutoff = (celltype_df %>% dplyr::filter(DNA_dam_class == 'damaged') %>% 
            pull(DNA_dam_val) %>% min() + 
          celltype_df %>% dplyr::filter(DNA_dam_class == 'undamaged') %>% 
            pull(DNA_dam_val) %>% max())/2

## add to Seurat object meta.data
obj = obj[,rownames(celltype_df)]
obj$DNA_dam_val = celltype_df$DNA_dam_val
obj$DNA_dam_class = celltype_df$DNA_dam_class


########################
# 3) average per sample
dam_per_sample = celltype_df%>% 
  mutate(Case = ss(ID, '-', 2)) %>% group_by(Case) %>% 
  mutate(numCell = n(), DNA_dam_prop = sum(DNA_dam_class =='damaged') / numCell) %>% 
  mutate_if(is.numeric, mean) %>% 
  ungroup() %>% distinct(Case, .keep_all = T) 

## stats
modBySample = lm(DNA_dam_val ~ DSM.IV.OUD+ Age + Sex + PMI + RIN + numCell,
                 data = dam_per_sample)
summary(modBySample) 

## no difference in proportion of Damaged cells in general
modBySample2 = glm(DNA_dam_prop ~ DSM.IV.OUD+ Age + Sex + PMI + RIN + numCell,
                 data = dam_per_sample, family = 'binomial')
summary(modBySample2) 

## plots
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamVal.perSample.pdf'), width = .75, height = 1)
ggplot(dam_per_sample, aes(x =DSM.IV.OUD, y = DNA_dam_val)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(pch = 21, alpha = 0.9, aes(x =DSM.IV.OUD, y = DNA_dam_val, fill = Sex)) + 
  geom_line(alpha = 0.9, aes(x =DSM.IV.OUD, y = DNA_dam_val, group=Pair, color = Sex)) + 
  scale_color_manual(values =c('pink', 'green')) +
  scale_fill_manual(values =c('pink', 'green')) +
  theme_bw(base_size = 5) +
  ylab('DNA Damage Score') + xlab('')+
  theme(legend.position = 'top', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()


##############################################
# 4) compute group average DNA damage ratios

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
         DNA_dam_prop = sum(DNA_dam_class =='damaged') / numCell) %>% 
  ungroup() %>% distinct(Case, Region, celltype4, .keep_all = T) %>% 
  ## z-normalize to rescale these numeric values for regression
  mutate_at(all_of(c('Age', 'PMI', 'RIN', 'numCell')), 
            ~ (. - mean(.))/sd(.)) %>% 
  # arrange(!grepl('D1', celltype4), !grepl('D2', celltype4), 
  #         !grepl('Int', celltype4), celltype4) %>% 
  group_by(celltype4) %>% mutate(tmp = mean(DNA_dam_val)) %>% 
  ungroup() %>% arrange(desc(DNA_dam_val)) %>%   
  mutate(celltype4 = factor(celltype4, unique(celltype4)))


## take a look
dam_per_cellxSample %>% count(Region, Case)
dam_per_cellxSample$DNA_dam_prop %>% summary()

## stats using mixed effect averages of cell type DNA dam scores per celltype, region, and patient
lmer(DNA_dam_val ~ celltype4 + celltype4:DSM.IV.OUD+ Age + Sex + 
       PMI + RIN + numCell + Region + (1|Case), 
     data = dam_per_cellxSample) %>% summary()

# Fixed effects:
#   Estimate Std. Error         df t value
# (Intercept)                         5.485e-02  1.683e-03  2.601e+01  32.600
# celltype4Endothelial               -1.014e-02  1.801e-03  1.055e+02  -5.626
# celltype4Microglia                 -1.458e-02  1.812e-03  1.056e+02  -8.046
# celltype4Astrocytes                -2.944e-02  1.817e-03  1.056e+02 -16.207
# celltype4Oligos                    -2.647e-02  2.402e-03  1.060e+02 -11.021
# celltype4Oligos_Pre                -3.526e-02  1.811e-03  1.056e+02 -19.472
# Age                                 1.127e-03  5.938e-04  7.318e+00   1.898
# SexM                               -1.758e-03  1.647e-03  8.432e+00  -1.067
# PMI                                -1.063e-03  8.280e-04  7.986e+00  -1.284
# RIN                                 5.600e-04  5.825e-04  6.325e+00   0.961
# numCell                             5.141e-04  5.605e-04  1.066e+02   0.917
# RegionPutamen                       5.195e-04  7.642e-04  1.082e+02   0.680
# celltype4Mural:DSM.IV.OUDOUD       -2.173e-04  2.027e-03  4.794e+01  -0.107
# celltype4Endothelial:DSM.IV.OUDOUD  5.448e-03  1.985e-03  4.571e+01   2.744
# celltype4Microglia:DSM.IV.OUDOUD    5.430e-03  1.948e-03  4.344e+01   2.788
# celltype4Astrocytes:DSM.IV.OUDOUD   2.900e-03  1.947e-03  4.339e+01   1.489
# celltype4Oligos:DSM.IV.OUDOUD       2.786e-03  1.945e-03  4.322e+01   1.432
# celltype4Oligos_Pre:DSM.IV.OUDOUD   4.340e-03  1.947e-03  4.338e+01   2.229
# Pr(>|t|)    
# (Intercept)                         < 2e-16 ***
#   celltype4Endothelial               1.52e-07 ***
#   celltype4Microglia                 1.35e-12 ***
#   celltype4Astrocytes                 < 2e-16 ***
#   celltype4Oligos                     < 2e-16 ***
#   celltype4Oligos_Pre                 < 2e-16 ***
#   Age                                 0.09768 .  
  # SexM                                0.31536    
  # PMI                                 0.23513    
  # RIN                                 0.37165    
  # numCell                             0.36109    
  # RegionPutamen                       0.49808    
  # celltype4Mural:DSM.IV.OUDOUD        0.91506    
  # celltype4Endothelial:DSM.IV.OUDOUD  0.00864 ** 
  # celltype4Microglia:DSM.IV.OUDOUD    0.00785 ** 
  # celltype4Astrocytes:DSM.IV.OUDOUD   0.14369    
  # celltype4Oligos:DSM.IV.OUDOUD       0.15925    
  # celltype4Oligos_Pre:DSM.IV.OUDOUD   0.03107 *  

modBySampleAndCell = lmer(DNA_dam_val ~ celltype4 + celltype4:DSM.IV.OUD+ Age + Sex + 
                            PMI + RIN + numCell + Region + (1|Case), 
                          data = dam_per_cellxSample) %>% tidy() %>% 
  filter(grepl('^celltype4', term) & grepl('OUD$', term)) %>% arrange(p.value) %>% 
  dplyr::select(-c('effect', 'group')) %>% 
  mutate(celltype = ss(term, '4|:', 2))
modBySampleAndCell %>% as.data.frame() %>% relocate('p.value', .after = 'term')

## plots
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamVal.perSampleByCell.pdf'), width = 4.5, height = 1)
ggplot(dam_per_cellxSample, 
       aes(x =DSM.IV.OUD, y = DNA_dam_val, fill = DSM.IV.OUD)) +
  geom_violin(size = .25) + geom_boxplot(width = 0.4, size = .25, outlier.shape = NA) + 
  # geom_point(pch = 21, alpha = 0.5, position = position_jitter(width = 0.3)) + 
  scale_fill_manual(values =c('white', 'red')) +
  scale_color_manual(values =c('white', 'red')) +
  facet_wrap(~celltype4, scale = 'fixed', nrow = 1)+
  theme_bw(base_size = 5) +
  ylab('DNA Damage Score') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()


###########################
# 5) export stats to table 
list('DNA_dam_score_by_sample' = summary(modBySample) %>% tidy(), 
     'DNA_dam_score_by_sampleCelltype' = modBySampleAndCell) %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_dnaDamVal.perSampleByCell.xlsx'))


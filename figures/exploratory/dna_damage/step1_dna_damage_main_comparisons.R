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
           "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
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
  
save_file = here(PLOTDIR, 'rdas', 'BU_OUD_Striatum_refined_all_SeuratObj_N22.meta.DNAdam.rds')
celltype_df %>% saveRDS(save_file)

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
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_dnaDamVal.perSample.pdf'), width = 1.25, height = 2)
ggplot(dam_per_sample, aes(x =DSM.IV.OUD, y = DNA_dam_val)) +
  geom_boxplot() + 
  geom_point(pch = 21, aes(x =DSM.IV.OUD, y = DNA_dam_val, fill = Sex)) + 
  geom_line(aes(x =DSM.IV.OUD, y = DNA_dam_val, group=Pair, color = Sex)) + 
  scale_color_manual(values =c('red', 'black')) +
  scale_fill_manual(values =c('red', 'black')) +
  theme_bw(base_size = 7) +
  ylab('DNA Damage Score') + xlab('Diagnosis')+
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
  arrange(!grepl('D1', celltype4), !grepl('D2', celltype4), 
          !grepl('Int', celltype4), celltype4) %>% 
  mutate(celltype4 = factor(celltype4, unique(celltype4)))

## take a look
dam_per_cellxSample %>% count(Region, Case)
dam_per_cellxSample$DNA_dam_prop %>% summary()

## stats using mixed effect averages of cell type DNA dam scores per celltype, region, and patient
modBySampleAndCell = lmer(DNA_dam_val ~ celltype4*DSM.IV.OUD+ Age + Sex + 
                            PMI + RIN + numCell + Region + (1|Case), 
                          data = dam_per_cellxSample) %>% tidy() %>% 
  filter(grepl('^celltype4', term) & grepl('OUD$', term)) %>% arrange(p.value) %>% 
  dplyr::select(-c('effect', 'group'))
modBySampleAndCell %>% as.data.frame() %>% relocate('p.value', .after = 'term')

## plots
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_dnaDamVal.perSampleByCell.pdf'), width = 4.5, height = 2)
ggplot(dam_per_cellxSample, aes(x =DSM.IV.OUD, y = DNA_dam_val, fill = DSM.IV.OUD)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(pch = 21, alpha = 0.5, position = position_jitter(width = 0.3)) + 
  scale_fill_manual(values =c('white', 'gray')) +
  scale_color_manual(values =c('white', 'black')) +
  facet_wrap(~celltype4, scale = 'free_y', nrow = 2)+
  theme_bw(base_size = 7) +
  ylab('DNA Damage Score') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()


###########################
# 5) export stats to table 
list('DNA_dam_score_by_sample' = summary(modBySample) %>% tidy(), 
     'DNA_dam_score_by_sampleCelltype' = modBySampleAndCell) %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'BU_OUD_Striatum_dnaDamVal.perSampleByCell.xlsx'))


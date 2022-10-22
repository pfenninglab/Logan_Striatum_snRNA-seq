## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(broom.mixed)
library(RColorBrewer)
library(rcartocolor)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## gene expression heatmap plots
library(AUCell)

## stats
library(tidymodels)
library(broom)
library(lme4)
library(lmerTest)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/AUCell_gene_set_activities'
PLOTDIR='figures/exploratory/AUCell_gene_set_activities/plots'
dir.create(here(PLOTDIR), showWarnings = F)



###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural/Fibroblast', 'Oligos', 'Oligos_Pre')
othertypes = c('Interneuron',  'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = c(carto_pal(length(othertypes) , "Vivid"))
names(othertypes_col) = othertypes

typecolors = c(subtypes_col, othertypes_col)
names(typecolors) = make.names(names(typecolors))

###################################################
## 1) read in the refined cell type metadata 
save_meta_fn = here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_SeuratObj_N22.metadata.rds')
meta = readRDS(file=save_meta_fn) %>% rownames_to_column('cellBarcode') 
table(meta$celltype3)


## load in the AUCell activities
save_fn = here(DATADIR, 'rdas', 'AUCell_Xue_OUDvsCONT_gainLoss_rhythmicity_refined_celltype_N22.rds')
cells_AUC = readRDS(save_fn)
cells_AUC_df =getAUC(cells_AUC) %>% t() %>% as.data.frame() %>% 
  rownames_to_column('cellBarcode') %>% 
  inner_join(x = meta, y = .)

cells_AUC_dflong = cells_AUC_df %>% 
  pivot_longer(cols = all_of(names(cells_AUC)), values_to = 'AUC', names_to = 'GeneSet') %>% 
  group_by(Case, Sex, Age, DSM.IV.OUD,  PMI, RIN, Region , celltype3, GeneSet) %>% 
  summarize(AUC= mean(AUC), numCell = n()) %>% ungroup()
  
###################################################
## 2) plot the differences in the AUCell gene sets
pdf(here(PLOTDIR, 'AUCell_Xue_rhythmicity_gainLoss_byCelltype3.NAc.pdf'), width = 6, height = 2.5)
p1 = ggplot(cells_AUC_dflong %>% filter(grepl("NAC", GeneSet)), 
            aes(x = DSM.IV.OUD, y = AUC)) + 
  geom_violin(size = .25, aes(fill = DSM.IV.OUD)) + 
  geom_boxplot(width = 0.4, size = .25, outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
  scale_fill_manual(values =c('white', '#00B57395'), guide = 'none') +
  # geom_text(data = modBySampleAndCell, size = 2, x = 1.5,
  #           aes(y = y_position + .02, label = annotation)) +
  facet_grid(GeneSet~celltype3, scale = 'free_x', space = 'fixed')+
  theme_bw(base_size = 4.5) +
  ylab('') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
p1
dev.off()


pdf(here(PLOTDIR, 'AUCell_Xue_rhythmicity_gainLoss_byCelltype3.DLPFC.pdf'), width = 6, height = 2.5)
p1 = ggplot(cells_AUC_dflong %>% filter(grepl("DLPFC", GeneSet)), 
            aes(x = DSM.IV.OUD, y = AUC)) + 
  geom_violin(size = .25, aes(fill = DSM.IV.OUD)) + 
  geom_boxplot(width = 0.4, size = .25, outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
  scale_fill_manual(values =c('white', '#00B57395'), guide = 'none') +
  # geom_text(data = modBySampleAndCell, size = 2, x = 1.5,
  #           aes(y = y_position + .02, label = annotation)) +
  facet_grid(GeneSet~celltype3, scale = 'free_x', space = 'fixed')+
  theme_bw(base_size = 4.5) +
  ylab('') + xlab('')+
  theme(legend.position = 'none', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), 
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
p1
dev.off()

##################################
## 3) mixed effect statistics 
stats_group = cells_AUC_dflong %>% 
  nest(data = -c(GeneSet)) %>% 
  mutate(
    test = map(data, ~ lmer(AUC ~ DSM.IV.OUD + celltype3 + + Age + Sex + 
                              PMI + RIN + numCell + Region + (1|Case), data = .x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  select(-data, -test) %>% 
  filter(grepl('DSM.IV.OUD', term)) %>% 
  mutate(FDR = p.adjust(p.value, 'fdr')) %>% 
  arrange(p.value) %>% 
  as.data.frame()

head(stats_group)









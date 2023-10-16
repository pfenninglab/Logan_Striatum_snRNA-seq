## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

## AUCell object for DNA damage estimates
library(AUCell)

## statistical analyses
library(tidymodels)
library(lme4)
library(broom.mixed)
library(lmerTest)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/compare_dna_damage_welch'
PLOTDIR='figures/exploratory/compare_dna_damage_welch'
dir.create(here(PLOTDIR, 'plots'), showWarnings = F)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F)
dir.create(here(PLOTDIR, 'rdas'), showWarnings = F)

########################################
# 1) read in Logan snRNA dataset to plot
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA') 
meta = obj_merged[[]] %>% rownames_to_column('cellBarcode') 

table(meta$celltype3)
table(meta$celltype4)

## load in the neuron-only AUCell scores
save_fn = here(DATADIR, 'rdas', 'AUCell_Welch_dna_dam_stages_refined_neurontype_N22.rds')
neuron_AUC = readRDS(save_fn)

## add the AUCell to the interneurons 
neuron_AUC_df =getAUC(neuron_AUC) %>% t() %>% as.data.frame() %>% 
  rownames_to_column('cellBarcode') %>% 
  inner_join(x = meta, y = .) %>% 
  filter(celltype4 == 'Interneurons') %>% 
  mutate(ID = factor(ID))

## linear modeling 
mod = lm(Stage1 ~ celltype3 + celltype3:DSM.IV.OUD+ Age + Sex + 
     PMI + RIN  + Region + (1|Case) , 
   data = neuron_AUC_df) %>% summary() %>% tidy() %>% 
  filter(grepl('celltype3Int', term), 
         grepl(':', term)) %>% 
  mutate(celltype3 = ss(term, '3|:', 2), 
         label = paste('p =', signif(p.value, 2)))

## plots
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamValStage1.byInterneurons.pdf'), width = 2.7, height = 1)
ggplot(neuron_AUC_df, aes(x =DSM.IV.OUD, y = Stage1)) +
  geom_violin(aes(fill = DSM.IV.OUD), size = .25, bw = .002) + 
  geom_boxplot(width = 0.4, size = .25, outlier.shape = NA) + 
  geom_jitter(aes(shape = Sex), size = .5, 
              position = position_jitterdodge(), alpha = .5)+
  geom_text(data = mod, aes(label = label, y = 0.045, x = 1), 
            size = 1.5) + 
  scale_fill_manual(values =c('white', '#C61B8C')) +
  scale_color_manual(values =c('white', '#C61B8C')) +
  scale_shape_manual(values =c(20, 18)) +
  facet_wrap(~celltype3, scale = 'fixed', nrow = 1)+
  theme_bw(base_size = 5) +
  ylab('NeuN+/gH2AX+ \nDamage Score') + 
  theme(legend.position = 'bottom', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), axis.title.x=element_blank(),
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-2,-10,-0,-10)) 
dev.off()





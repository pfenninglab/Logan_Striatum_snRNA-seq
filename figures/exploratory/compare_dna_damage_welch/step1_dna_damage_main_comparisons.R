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

#########################################
# 1) load in the DNA damage estimates

## load in the glia only AUCell scores
save_meta_fn = here('data/tidy_data/AUCell_gene_set_activities/rdas/BU_OUD_Striatum_refined_all_SeuratObj_N22.metadata.rds')
meta = readRDS(file=save_meta_fn) %>% rownames_to_column('cellBarcode') 
table(meta$celltype3)

save_fn = here(DATADIR, 'rdas', 'AUCell_Welch_dna_dam_stages_refined_gliatype_N22.rds')
glia_AUC = readRDS(save_fn)

## calculate adaptive thresholds to call a "damaged" glial cell
AUCell_thr_fn = here(PLOTDIR, 'plots', 'AUCell_welch_dna_damage_scores_thresholds.glia.pdf')
pdf(AUCell_thr_fn, onefile = T)
glia_assignment <- AUCell_exploreThresholds(glia_AUC, plotHist=T, assign=T, nCores = 1) 
dev.off()

lapply(glia_assignment, '[[', 'aucThr')

glia_AUC_df = getAUC(glia_AUC) %>% t() %>% as.data.frame() %>% 
  rownames_to_column('cellBarcode') %>% 
  inner_join(x = meta, y = .) %>% 
  filter(!grepl('^D', celltype3), celltype3!= 'Interneuron' ) %>% 
  mutate(celltype3 = droplevels(celltype3), 
         Stage2_class = cellBarcode %in% glia_assignment$Stage2$assignment)
table(glia_AUC_df$celltype3)
table(glia_AUC_df$Stage2_class)

## load in the neuron-only AUCell scores
save_fn = here(DATADIR, 'rdas', 'AUCell_Welch_dna_dam_stages_refined_neurontype_N22.rds')
neuron_AUC = readRDS(save_fn)

## calculate adaptive thresholds to call a "damaged" neuronal cell
AUCell_thr_fn2 = here(PLOTDIR, 'plots', 'AUCell_welch_dna_damage_scores_thresholds.neuron.pdf')
pdf(AUCell_thr_fn2, onefile = T)
neuron_assignment <- AUCell_exploreThresholds(neuron_AUC, plotHist=T, assign=T, nCores = 1) 
dev.off()

lapply(neuron_assignment, '[[', 'aucThr')

neuron_AUC_df =getAUC(neuron_AUC) %>% t() %>% as.data.frame() %>% 
  rownames_to_column('cellBarcode') %>% 
  inner_join(x = meta, y = .) %>% 
  mutate(celltype3 = droplevels(celltype3), 
         Stage1_class = cellBarcode %in% neuron_assignment$Stage1$assignment)
table(neuron_AUC_df$celltype3)
table(neuron_AUC_df$Stage1_class)

##############################################
# 2) average neuron scores per sample neurons
neuron_dam_per_sample = neuron_AUC_df %>% group_by(Case) %>% 
  mutate(Stage1= mean(Stage1), 
         numCell = n()) %>% ungroup() %>% 
  distinct(DSM.IV.OUD, Age, Sex, PMI, RIN, Region, numCell, .keep_all = T)

## stats
neuron_modBySample_s1 = lm(Stage1 ~ DSM.IV.OUD + Age + Sex + PMI + RIN + numCell,
                    data = neuron_dam_per_sample) 
summary(neuron_modBySample_s1) 

# lm(formula = Stage1 ~ DSM.IV.OUD + Age + Sex + PMI + RIN + numCell, 
#    data = neuron_dam_per_sample)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -2.469e-03 -7.856e-04 -4.803e-05  6.431e-04  1.667e-03 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)    1.333e-02  3.270e-03   4.078  0.00099 ***
#   DSM.IV.OUDOUD  1.303e-03  5.980e-04   2.179  0.04569 *  
#   Age           -1.266e-05  3.751e-05  -0.337  0.74048    
#   SexM          -1.374e-03  8.524e-04  -1.612  0.12771    
#   PMI            2.437e-04  1.222e-04   1.994  0.06470 .  
#   RIN           -1.558e-04  4.147e-04  -0.376  0.71241    
#   numCell       -2.067e-06  8.540e-07  -2.420  0.02869 *  

## plots
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum.neuron.dnaDamStage1.perSample.pdf'), width = .75, height = 1)
ggplot(neuron_dam_per_sample, aes(x =DSM.IV.OUD, y = Stage1, fill = DSM.IV.OUD)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = 0.9, aes(x =DSM.IV.OUD, y = Stage1, shape = Sex)) + 
  scale_fill_manual(values =c('white', '#C61B8C')) +
  scale_shape_manual(values =c(21, 22)) +
  theme_bw(base_size = 5) +
  ylab('NeuN+/gH2AX+\nDamage Score') +
  theme(legend.position = 'none', axis.title.x=element_blank()) 
dev.off()


pdf(here(PLOTDIR, 'plots', 'OUD_Striatum.neuron.dnaDamStage1.perSampleByPMI.pdf'), width = 2.25, height = 1)
ggplot(neuron_dam_per_sample, aes(x =PMI, y = Stage1, color = DSM.IV.OUD)) +
  geom_smooth(method="lm", se=F, fullrange=T, level=0.95)+
  geom_point(alpha = 0.9, aes(shape = Sex)) + 
  scale_shape_manual(values =c(21, 22)) +
  scale_color_manual(values =c('black', '#C61B8C')) +
  theme_bw(base_size = 5) + facet_wrap(~Sex, scales = 'free_x') + 
  ylab('NeuN+/gH2AX+\nDamage Score') + xlab("PMI") +
  theme(legend.position = 'none') 
dev.off()


####################################################
# 4) average glia DNA damage scores per sample glia
glia_dam_per_sample = glia_AUC_df %>% group_by(Case) %>% 
  mutate(Stage2 = mean(Stage2), 
         numCell = n()) %>% ungroup() %>% 
  distinct(DSM.IV.OUD, Age, Sex, PMI, RIN, Region, numCell, .keep_all = T)

glia_modBySample_s2 = lm(Stage2 ~ DSM.IV.OUD+ Age + Sex + PMI + RIN + numCell,
                    data = glia_dam_per_sample)
summary(glia_modBySample_s2) 

# lm(formula = Stage2 ~ DSM.IV.OUD + Age + Sex + PMI + RIN + numCell, 
#    data = glia_dam_per_sample)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0052424 -0.0016133  0.0004051  0.0017710  0.0032343 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)    1.511e-01  7.952e-03  19.000 6.62e-12 ***
#   DSM.IV.OUDOUD  3.531e-03  1.341e-03   2.633 0.018807 *  
#   Age            9.827e-05  1.447e-04   0.679 0.507321    
#   SexM           6.032e-03  2.128e-03   2.835 0.012543 *  
#   PMI           -1.347e-03  2.654e-04  -5.077 0.000136 ***
#   RIN            1.128e-03  1.170e-03   0.964 0.350191    
#   numCell        1.643e-07  4.158e-07   0.395 0.698261    


pdf(here(PLOTDIR, 'plots', 'OUD_Striatum.glia.dnaDamStage2.perSample.pdf'), width = .75, height = 1)
ggplot(neuron_dam_per_sample, aes(x =DSM.IV.OUD, y = Stage2, fill = DSM.IV.OUD)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(alpha = 0.9, aes(x = DSM.IV.OUD, y = Stage2, shape = Sex)) + 
  scale_shape_manual(values =c(21, 22)) +
  scale_fill_manual(values =c('white', '#2A62AD')) +
  theme_bw(base_size = 5) +
  ylab('NeuN-/gH2AX+\nDamage Score') +
  theme(legend.position = 'none', axis.title.x=element_blank()) 
dev.off()


pdf(here(PLOTDIR, 'plots', 'OUD_Striatum.glia.dnaDamStage2.perSampleByPMI.pdf'), width = 2.25, height = 1)
ggplot(neuron_dam_per_sample, aes(x =PMI, y = Stage2, color = DSM.IV.OUD)) +
  geom_smooth(method="lm", se=F, fullrange=T, level=0.95)+
  geom_point(alpha = 0.9, aes(shape = Sex)) + 
  scale_shape_manual(values =c(21, 22)) +
  scale_color_manual(values =c('black', '#2A62AD')) +
  theme_bw(base_size = 5) + facet_wrap(~Sex, scales = 'free_x') + 
  ylab('NeuN-/gH2AX+\nDamage Score') +xlab("PMI") +
  theme(legend.position = 'none') 
dev.off()


######################################################################
# 5) compute neuron type average DNA damage ratios for neuronal damage

## stratified by cell type and sample
dam_per_neuronxSample = neuron_AUC_df %>% 
  arrange(celltype3) %>% 
  group_by(Case, Region, celltype3) %>% 
  mutate_if(is.numeric, mean) %>% 
  mutate(numCell = n()) %>% 
  ungroup() %>% distinct(Case, Region, celltype3, .keep_all = T) %>% 
  ## z-normalize to rescale these numeric values for regression
  mutate_at(all_of(c('Age', 'PMI', 'RIN', 'numCell')), 
            ~ (. - mean(.))/sd(.)) %>% 
  group_by(celltype3) %>%   
  mutate(celltype3 = factor(celltype3, unique(celltype3)))

lm(Stage1 ~ celltype3 + celltype3:DSM.IV.OUD+ Age + Sex + 
       PMI + RIN + numCell + Region , 
     data = dam_per_neuronxSample) %>% summary()

## stage 1 in glia
modBySampleAndcellStage1 = lm(Stage1 ~ celltype3 + celltype3:DSM.IV.OUD+ Age + Sex + 
                            PMI + RIN + numCell + Region + (1|Case), 
                          data = dam_per_neuronxSample) %>% tidy() %>% 
  filter(grepl('^celltype3', term) & grepl('OUD$', term)) %>% arrange(p.value) %>% 
  mutate(celltype = ss(term, '3|:', 2))

# lm(formula = Stage1 ~ celltype3 + celltype3:DSM.IV.OUD + Age + 
#      Sex + PMI + RIN + numCell + Region, data = dam_per_neuronxSample)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0037683 -0.0012404 -0.0003683  0.0007216  0.0085088 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept)                          0.0126353  0.0008627  14.646  < 2e-16 ***
#   celltype3D2.Matrix                   0.0006177  0.0010243   0.603  0.54772    
#   celltype3D1.Striosome               -0.0004952  0.0010913  -0.454  0.65088    
#   celltype3D2.Striosome                0.0004463  0.0010868   0.411  0.68211    
#   celltype3D1.D2.Hybrid                0.0011699  0.0011144   1.050  0.29606    
#   celltype3Interneuron                 0.0033980  0.0010913   3.114  0.00234 ** 
#   Age                                  0.0006139  0.0002172   2.826  0.00556 ** 
#   SexM                                -0.0016553  0.0006430  -2.574  0.01132 *  
#   PMI                                  0.0005207  0.0003208   1.623  0.10726    
#   RIN                                 -0.0004990  0.0002094  -2.383  0.01882 *  
#   numCell                             -0.0004463  0.0003520  -1.268  0.20737    
#   RegionPutamen                        0.0003189  0.0004686   0.681  0.49749    
#   celltype3D1.Matrix:DSM.IV.OUDOUD     0.0013839  0.0009866   1.403  0.16342    
#   celltype3D2.Matrix:DSM.IV.OUDOUD     0.0007931  0.0009911   0.800  0.42526    
#   celltype3D1.Striosome:DSM.IV.OUDOUD  0.0009064  0.0009816   0.923  0.35778    
#   celltype3D2.Striosome:DSM.IV.OUDOUD -0.0001250  0.0009817  -0.127  0.89888    
#   celltype3D1.D2.Hybrid:DSM.IV.OUDOUD  0.0001878  0.0009818   0.191  0.84862    
#   celltype3Interneuron:DSM.IV.OUDOUD   0.0021008  0.0009816   2.140  0.03448 * 


## plots
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamValStage1.perSampleByCell.pdf'), width = 2.7, height = 1)
ggplot(dam_per_neuronxSample, 
       aes(x =DSM.IV.OUD, y = Stage1, fill = DSM.IV.OUD)) +
  geom_violin(size = .25, bw = .002) + 
  geom_boxplot(width = 0.4, size = .25, outlier.shape = NA) + 
  geom_jitter(aes(shape = Sex), size = .5, 
              position = position_jitterdodge(), alpha = .5)+
  scale_fill_manual(values =c('white', '#C61B8C')) +
  scale_color_manual(values =c('white', '#C61B8C')) +
  scale_shape_manual(values =c(21, 22)) +
  facet_wrap(~celltype3, scale = 'fixed', nrow = 1)+
  theme_bw(base_size = 5) +
  ylab('NeuN+/gH2AX+ \nDamage Score') + 
  theme(legend.position = 'bottom', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), axis.title.x=element_blank(),
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-2,-10,-0,-10)) 
dev.off()



################################################################
# 6) compute cell type average DNA damage ratios for glial damage

## stratified by cell type and sample
dam_per_gliaxSample = glia_AUC_df %>% 
  arrange(celltype3) %>% 
  group_by(Case, Region, celltype3) %>% 
  mutate_if(is.numeric, mean) %>% 
  mutate(numCell = n()) %>% 
  ungroup() %>% distinct(Case, Region, celltype3, .keep_all = T) %>% 
  ## z-normalize to rescale these numeric values for regression
  mutate_at(all_of(c('Age', 'PMI', 'RIN', 'numCell')), 
            ~ (. - mean(.))/sd(.)) %>% 
  group_by(celltype3) %>%   
  mutate(celltype3 = factor(celltype3, unique(celltype3)))

## take a look
dam_per_gliaxSample %>% count(Region, Case)

## stats using mixed effect averages of cell type DNA dam scores per celltype, region, and patient
lm(Stage2 ~ celltype3 + celltype3:DSM.IV.OUD+ Age + Sex + 
     PMI + RIN + numCell + Region , 
   data = dam_per_gliaxSample) %>% summary()

# lm(formula = Stage2 ~ celltype3 + celltype3:DSM.IV.OUD + Age + 
#      Sex + PMI + RIN + numCell + Region, data = dam_per_gliaxSample)
# 
# Residuals:
#   Min        1Q    Median        3Q       Max 
# -0.037900 -0.004827 -0.000574  0.004778  0.040811 
# 
# Coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
#   (Intercept)                         0.1452904  0.0040240  36.106  < 2e-16 ***
#   celltype3Endothelial                0.0535110  0.0051669  10.356  < 2e-16 ***
#   celltype3Microglia                  0.0332393  0.0051282   6.482 2.46e-09 ***
#   celltype3Mural                      0.0329734  0.0051707   6.377 4.08e-09 ***
#   celltype3Oligos                    -0.0020299  0.0064234  -0.316  0.75257    
#   celltype3Oligos_Pre                -0.0441612  0.0051289  -8.610 4.92e-14 ***
#   Age                                 0.0036937  0.0011169   3.307  0.00126 ** 
#   SexM                               -0.0030465  0.0032211  -0.946  0.34627    
#   PMI                                -0.0042015  0.0016052  -2.617  0.01007 *  
#   RIN                                -0.0005167  0.0010492  -0.492  0.62334    
#   numCell                             0.0022630  0.0019921   1.136  0.25836    
#   RegionPutamen                       0.0056320  0.0021211   2.655  0.00907 ** 
#   celltype3Astrocytes:DSM.IV.OUDOUD   0.0063944  0.0049343   1.296  0.19765    
#   celltype3Endothelial:DSM.IV.OUDOUD  0.0101785  0.0049300   2.065  0.04125 *  
#   celltype3Microglia:DSM.IV.OUDOUD    0.0090970  0.0049365   1.843  0.06798 .  
#   celltype3Mural:DSM.IV.OUDOUD       -0.0016355  0.0050338  -0.325  0.74586    
#   celltype3Oligos:DSM.IV.OUDOUD       0.0054766  0.0049284   1.111  0.26882    
#   celltype3Oligos_Pre:DSM.IV.OUDOUD   0.0080470  0.0049337   1.631  0.10567    

## stage 2 in glia
modBySampleAndcellStage2 = lm(Stage2 ~ celltype3 + celltype3:DSM.IV.OUD+ Age + Sex + 
                                PMI + RIN + numCell + Region + (1|Case), 
                              data = dam_per_gliaxSample) %>% tidy() %>% 
  filter(grepl('^celltype3', term) & grepl('OUD$', term)) %>% arrange(p.value) %>% 
  mutate(celltype = ss(term, '3|:', 2))


## plots
pdf(here(PLOTDIR, 'plots', 'OUD_Striatum_dnaDamValStage2.perSampleByCell.pdf'), width = 2.7, height = 1)
ggplot(dam_per_gliaxSample, aes(x =DSM.IV.OUD, y = Stage1, fill = DSM.IV.OUD)) +
  geom_violin(size = .25, bw = .002) + 
  geom_boxplot(width = 0.4, size = .25, outlier.shape = NA) + 
  geom_jitter(aes(shape = Sex), size = .5, 
              position = position_jitterdodge(), alpha = .5)+
  scale_fill_manual(values =c('white', '#2A62AD')) +
  scale_color_manual(values =c('white', '#2A62AD')) +
  scale_shape_manual(values =c(21, 22)) +
  facet_wrap(~celltype3, scale = 'fixed', nrow = 1)+
  theme_bw(base_size = 5) +
  ylab('NeuN-/gH2AX+\nDamage Score') +
  theme(legend.position = 'bottom', legend.title = element_blank(),
        legend.spacing.y = unit(.5, 'cm'), axis.title.x=element_blank(),
        legend.key.size = unit(.2, "cm"), 
        legend.box.margin=margin(-5,-10,-5,-10)) 
dev.off()

#############################
# 7) export to spreadsheet
dfList = list(
  'Stage1_NeuN+_gH2AX+_neuron' = 
    modBySampleAndcellStage1 %>% as.data.frame() %>% relocate('celltype', .after = 'term'), 
  'Stage2_NeuN-_gH2AX+_glia' = 
    modBySampleAndcellStage2 %>% as.data.frame() %>% relocate('celltype', .after = 'term'))

stats_fn = here(PLOTDIR, 'tables', 'OUD_Striatum_dnaDam_statistics.perSampleByCell.xlsx')
writexl::write_xlsx(dfList, stats_fn)


ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggsci)
library(here)

DATADIR='data/tidy_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure5_sex_interaction_degs'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

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

in2mm<-25.4
my_theme = theme_classic(base_size = 6)

#################################################
# 1) load in the DEG table from the big analyses
rdasDir =file.path(DATADIR,'differential_expression_analysis', 'rdas')

## DEGs of OUD vs. Control in either females or male
res_OUDwinSex = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds') %>% readRDS() 

######################################
## 2) count DEGs within each grouping
alpha = 0.05
df = res_OUDwinSex %>% lapply(function(x){
  ## calculate a number that gives both statistical signif + effect size
  x = x %>% mutate(tmpF = -log10(adj.P.Val.Between_SexF) * sign(logFC_SexF),
         tmpM = -log10(adj.P.Val.Between_SexM) * sign(logFC_SexM),
         ## interaction_score is based on a mean-difference MA-plot
         ## this will score genes (most different vs. most similar of OUD DEGs in M or F)
         interaction_score = (tmpF - tmpM)/ abs(tmpM + tmpF)) %>% 
    filter(adj.P.Val.Between_SexF < alpha | adj.P.Val.Between_SexM < alpha)
    
  out = data.frame("Female.Unique" = sum(x$interaction_score > 0,  na.rm = T),
                   "Male.Unique" = sum(x$interaction_score < 0,  na.rm = T))
  return(out)
}) %>% data.table::rbindlist(idcol = 'celltype')

df_long = df %>% pivot_longer(cols = c('Male.Unique', 'Female.Unique'), names_to = 'Direction', 
                              values_to = 'numDEG') %>% 
  mutate(numDEG = ifelse(Direction =='Female.Unique', -numDEG, numDEG),
         Direction = factor(Direction, c('Female.Unique', 'Male.Unique')), 
         celltype_class = case_when(
           grepl('All|Neuron|Glia', celltype) ~ 'Class',
           grepl('^D|^Int', celltype)~ 'Neuron', 
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Class', 'Neuron', 'Glia')), 
         celltype = factor(celltype,  names(typecolors)))


##########################
## 2) make the miami plot
plot_fn = here(PLOTDIR,'plots', 'figure5_mumDEG_sexInteraction_miamiPlot.pdf')
pdf(plot_fn, height = 80/in2mm, width = 50/in2mm)
ggplot(df_long, aes(x = celltype, y = numDEG, fill = celltype))+
  geom_bar(stat = 'identity')+
  geom_text(size = 1.5, aes(label=abs(numDEG)), hjust="inward") + 
  facet_grid( celltype_class ~ Direction, scales = 'free', space = 'free')+
  scale_fill_manual(values = typecolors) +
  scale_x_discrete(limits=rev) +
  scale_y_continuous(labels = abs) + # so negative sign doesn't show
  my_theme + coord_flip() + 
  ylab("# of DEGs, FDR < 0.05") + 
  theme(legend.position = 'none', 
        axis.title.y = element_blank(), 
        axis.text.x = element_blank())
dev.off()


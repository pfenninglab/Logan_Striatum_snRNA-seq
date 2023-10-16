ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggsci)
library(ggh4x)
library(here)

DATADIR='data/tidy_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure2_the_glia_degs_and_cell_states'
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
alpha = 0.05

#################################################
# 1) load in the DEG table from the big analyses
rdasDir =file.path(DATADIR,'differential_expression_analysis', 'rdas')
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds')
res = readRDS(save_res_fn)

# number of cell types
length(res) 

# number of unique genes that are DEGs across any cell type
res %>% lapply(function(x) x %>% filter(adj.P.Val.Between < alpha)) %>% 
  rbindlist() %>% filter(!duplicated(gene)) %>% nrow()


######################################
## 2) count DEGs within each grouping
df = res %>% lapply(function(x){
  out = data.frame(numDiff = sum(x$adj.P.Val.Between < alpha, na.rm = T), 
                   "Up.Regulated" = sum(x$adj.P.Val.Between < alpha & x$logFC > 0,  na.rm = T),
                   "Down.Regulated" = sum(x$adj.P.Val.Between < alpha & x$logFC < 0,  na.rm = T))
  return(out)
}) %>% data.table::rbindlist(idcol = 'celltype')

## mean and SE of DEGs per cell type
df %>% pull(numDiff) %>% mean()
df %>% pull(numDiff) %>% sd() /sqrt(length(res))

df_long = df %>% pivot_longer(cols = c('Up.Regulated', 'Down.Regulated'), names_to = 'Direction', 
                              values_to = 'numDEG') %>% 
  mutate(numDEG = ifelse(Direction =='Down.Regulated', -numDEG, numDEG),
         Direction = factor(Direction, c( 'Up.Regulated', 'Down.Regulated')), 
         celltype_class = case_when(
           grepl('All|Neuron|Glia', celltype) ~ 'Class',
           grepl('^D|^Int', celltype)~ 'Neuron', 
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Class', 'Neuron', 'Glia')), 
         celltype = factor(celltype,  names(typecolors))) %>% 
  dplyr::filter(celltype_class!= 'Class')


##########################
## 2) make the miami plot
plot_fn = here(PLOTDIR,'plots', 'figure2_mumDEG_celltypes_miamiPlot_glia3.pdf')
pdf(plot_fn, height = 45/in2mm, width = 80/in2mm)
ggplot(df_long, aes(x = celltype, y = numDEG, fill = celltype))+
  geom_bar(stat = 'identity')+
  geom_text(size = 2, aes(label=abs(numDEG)),vjust="inward") + 
  facet_nested( Direction ~  celltype_class, scales = 'free', space = 'free')+
  scale_fill_manual(values = typecolors) + my_theme + 
  scale_y_continuous(labels = abs) + # so negative sign doesn't show
  ylab("# of DEGs, FDR < 0.05") + 
  theme(legend.position = 'none', axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 15, vjust = 1, hjust=1))
dev.off()


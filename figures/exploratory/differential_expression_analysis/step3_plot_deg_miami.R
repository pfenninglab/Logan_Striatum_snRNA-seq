## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)
library(svglite)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'
h5Dir =here(DATADIR, 'HDF5Array')
dir.create(h5Dir, showWarnings = F)

PLOTDIR='figures/exploratory/differential_expression_analysis/plots'
dir.create(PLOTDIR, showWarnings = F)


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


## load in the DEG table
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'BU_OUD_Striatum_voom_limmaRes_bigModel_N22.rds')
res = readRDS(save_res_fn)

alpha = 0.05
df = res %>% lapply(function(x){
  out = data.frame(numDiff = sum(x$p_adj < alpha, na.rm = T), 
                   UpReg = sum(x$p_adj < alpha & x$logFC > 0,  na.rm = T),
                   DnReg = sum(x$p_adj < alpha & x$logFC < 0,  na.rm = T))
  return(out)
}) %>% rbindlist(idcol = 'celltype') %>% 
  filter(!celltype %in% c( "DSM.IV.OUD", "RegionPutamen", "SexM" ) )

df_long = df %>% pivot_longer(cols = c('UpReg', 'DnReg'), names_to = 'Direction', 
                              values_to = 'numDEG') %>% 
  mutate(numDEG = ifelse(Direction =='DnReg', -numDEG, numDEG),
         Direction = factor(Direction, c('UpReg', 'DnReg')), 
         celltype_class = case_when(grepl('^D|^In', celltype)~ 'Neuron', 
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia')), 
         celltype = factor(celltype, names(typecolors)))

## make the plot
plot_fn = here(PLOTDIR, 'BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.miamiNumDEG.pdf')
pdf(plot_fn, height = 4, width = 8)
ggplot(df_long, aes(x = celltype, y = numDEG, fill = celltype))+
  geom_bar(stat = 'identity')+
  geom_text(aes(label=abs(numDEG)), vjust="inward") + 
  facet_grid(Direction ~ celltype_class, scales = 'free', space = 'free')+
  scale_fill_manual(values = typecolors) +
  scale_y_continuous(labels = abs) + # so negative sign doesn't show
  theme_bw(base_size = 10) + 
  theme(axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))
dev.off()

##alt save
plot_fn = here(PLOTDIR, 'BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.miamiNumDEG.svg')
ggsave(file = "BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.miamiNumDEG.svg", plot = image, width = 15, height = 4)


plot_fn = here(PLOTDIR, 'BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.miamiNumDEG.svg')
svglite(
  filename = "BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.miamiNumDEG.svg",
  width = 15,
  height = 4,
  bg = "white",
  pointsize = 12,
  standalone = TRUE,
  system_fonts = list(),
  user_fonts = list(),
  web_fonts = list(),
  id = NULL,
  fix_text_size = TRUE,
  scaling = 1,
  always_valid = FALSE,
  file
)

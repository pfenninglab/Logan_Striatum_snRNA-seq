ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(SingleCellExperiment)
library(ComplexHeatmap)
library(rcartocolor)
library(data.table)
library(tidyverse)
library(ggsci)
library(here)
library(readxl)
library(ggh4x)
library(wesanderson)
library(cividis) #  devtools::install_github("marcosci/cividis")


DATADIR='data/tidy_data/differential_expression_analysis'

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

dx_col = setNames(pal_npg(palette = c("nrc"), alpha = 1)(2), c('CTL', 'OUD'))

in2mm<-25.4
my_theme = theme_classic(base_size = 5)

#############################################
# 1) load in the pseudobulk of TF regulators
alpha = 0.05

## read in the saved pseudo bulk TF regulon output 
stats_group  = here('figures/exploratory/AUCell_gene_set_activities', 'tables', 'OUD_Striatum_pyscenic_differential_TF_modules.byCelltype.xlsx') %>% 
  readxl::read_xlsx() %>% 
  filter(celltype3 %in% c('Astrocytes','Endothelial', 'Microglia', 'Mural', 'Oligos', 'Oligos_Pre')) %>% 
  mutate(isSignif = p.adj < alpha )

signif_TF = stats_group %>% filter(isSignif ) %>% pull(TF) %>% unique()
signif_cell = stats_group %>% filter(isSignif ) %>% pull(celltype3) %>% unique()

## reorder the TFs by the difference in OUD
stats_group2 = stats_group %>% filter(TF %in% signif_TF, celltype3 %in% signif_cell) %>% 
  group_by(TF) %>% mutate(tmp = mean(statistic)) %>% ungroup() %>% 
  arrange(tmp) %>% mutate(TF = factor(TF, unique(TF)))

stats_group2 %>% dplyr::select(-tmp) %>%
  writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_pyscenic_differential_TF_modules.byCelltype.sourceData.xlsx'))

#############################################
# 2) load in the pseudobulk of TF regulators
plot_tf_enrichment_fn = here(PLOTDIR, 'plots', 'figure2_TF_regulon_glia.heatmap.pdf')
pdf(plot_tf_enrichment_fn, height = 70/in2mm, width = 59/in2mm)
ggplot(stats_group2, aes(x = celltype3, y = TF, fill = statistic)) + 
  geom_tile(aes(color = isSignif), linewidth = .5, width = .95, height = .9) + 
  scale_fill_distiller(palette = "PiYG",direction =1, 't-statistic')+
  scale_color_manual(values = c('#00000000', 'black'), 'FDR < 0.05')+
  guides(color = guide_legend(override.aes = list(fill = "white"))) +
  my_theme + xlab('Cell type') + 
  ylab('Differential TF-gene regulatory module in OUD')
dev.off()


ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(ComplexHeatmap)
library(rcartocolor)
library(ggsci)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(SingleCellExperiment)
library(future)
library(here)

## function for pseudobulk based findAllMarkers
library(CVRCFunc) # devtools::install_github('mgildea87/CVRCFunc')
library(swfdr)

DATADIR='data/tidy_data/differential_expression_analysis'

## make for this subdirs
PLOTDIR='figures/explanatory/figure1_the_experiment'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

#######################################################
# 0) Seurat uses the future package for parallelization
plan("sequential")
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

in2mm<-25.4
alpha = 0.05
my_theme = theme_classic(base_size = 5)

subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes %>% make.names()

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'All', 'Neuron', 'Glia',
               'Astrocytes', 'Endothelial', 'Microglia', 'Mural', 'Oligos', 
               'Oligos_Pre', 'Interneurons')
othertypes_col = c(carto_pal(4, "Safe"), 
                   mypal = pal_npg('nrc')(3),
                   carto_pal(length(othertypes) -7 , "Vivid"))
names(othertypes_col) = othertypes
othertypes_col = othertypes_col[-c(1:4)]
typecolors = c(subtypes_col, othertypes_col[c(10, 4:9)])

########################################
# 1) read in Logan snRNA dataset to plot
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA') 
names(obj_merged[[]] )
obj_merged$celltype4 = ifelse(grepl('Int', obj_merged$celltype3), 
                              'Interneurons', obj_merged$celltype3) %>% make.names()
table(obj_merged$celltype4 )

cells = c( 'D1.Matrix', 'D1.Striosome', 'D2.Matrix', 'D2.Striosome','D1.D2.Hybrid')
obj_neuron = subset(obj_merged, celltype4 %in% cells)

########################################
# 2) extract marker genes 
genes = c('DRD1', 'DRD2', 'FOXP2', 'OTOF', 'CASZ1', 'CACNG5')

# extract these genes and plot co-expression
df_expr = cbind(obj_neuron[[]], FetchData(obj_neuron, genes)) %>% 
  mutate(celltype4 = factor(celltype4, cells))

df_expr_long = df_expr %>% 
  pivot_longer(cols = all_of(genes), names_to = 'gene') %>% 
  mutate(gene = factor(gene, genes))

plot_fn = here(PLOTDIR, 'plots', 'figure1_DRD1_DRD2_scatterplot.pdf')
pdf(plot_fn, height = 11, width = 8.5, onefile = T)

ggplot(df_expr, aes(x = DRD1, y = DRD2)) + 
  geom_hex() +
  geom_abline(slope = 1, color = 'red') + 
  geom_hline(yintercept = 5) + geom_vline(xintercept = 2 ) +
  facet_wrap(~celltype4, ncol = 2) +
  theme(legend.position = 'none') + 
  theme_bw() + 
  viridis::scale_color_viridis() +
  scale_x_continuous( limits = c(NA, 90)) 
dev.off()



plot_fn = here(PLOTDIR, 'plots', 'figure1_FOXP2_CASZ1_OTOF_D1_D2_violin_plot.pdf')
pdf(plot_fn, height = 11, width = 8.5, onefile = T)
ggplot(df_expr_long, aes(x = celltype4, y = value)) + 
  geom_violin(aes(fill = celltype4)) + 
  facet_wrap(~gene, ncol = 2, scales = 'free_x') + theme_bw() + 
  scale_fill_manual(values = subtypes_col) + 
  scale_y_continuous( trans = 'log10') +
  coord_flip() + 
  theme(legend.position = 'none') 
dev.off()

########################################

## make marker gene plot
plot_fn = here(PLOTDIR, 'plots', 'figure1_pseudobulk_celltype_markerGene_heatmap.pdf')
pdf(plot_fn, height = 60/in2mm, width = 90/in2mm, onefile = T)
ggplot(df, aes(x = groupings, y = gene, fill = value))+
  geom_tile() + 
  scale_fill_distiller(palette = "Spectral", 'Scaled Expression') + 
  scale_y_discrete(limits=rev) + 
  facet_grid(~celltype, scales = 'free_x') +
  my_theme + xlab('Biological Sample') +
  theme_classic(base_size = 4) +
  theme( legend.key.height = unit(2, 'mm'), legend.position = 'bottom', 
         axis.text.x = element_blank(), axis.ticks
         = element_blank(), 
         strip.text.x = element_text(size = 2.45),
         legend.box.margin=margin(-2,-2,-2,-2))
dev.off()



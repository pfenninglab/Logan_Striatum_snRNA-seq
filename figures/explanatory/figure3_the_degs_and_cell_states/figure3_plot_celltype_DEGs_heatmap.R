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
PLOTDIR='figures/explanatory/figure3_the_degs_and_cell_states'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)
source(here('code','final_code','Rutils', 'plot_logCPM_pseudobulk.R'))

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
my_theme = theme_classic(base_size = 6)

###################################
# 1) expression matrix to regress
v = readRDS(here(DATADIR, 'rdas', 'voomLimma_norm_object_N222pb.rds'))
designSV = readRDS(here(DATADIR, 'rdas', 'bigModelFitSVA_designMatrix.rds'))

## grab the phenotype data
df = here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds') %>% 
  readRDS() %>% colData() %>% as.data.frame() %>% 
  dplyr::select(ID:CaseIdx) %>% 
  mutate(celltype3 = celltype3 %>% make.names %>% factor(names(typecolors)),
         celltype3 = droplevels(celltype3), 
         celltype_class = case_when(grepl('^D|^Int', celltype3)~ 'Neuron',  
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia')))
df = df[colnames(v), ]

## move the cell type by diagnoses before all the covariables
designSV = designSV %>% as.data.frame() %>% 
  dplyr::relocate(contains('celltype_dx'), .before = everything()) %>% as.matrix()
num_regress = sum(grepl('celltype_dx', colnames(designSV)))/3*2

## regress everything except the cell type by diagnoses & interaction effects
## the $E slot is log2 CPM, https://rdrr.io/bioc/limma/man/EList.html
# y_clean = 2 ^ cleaningY(v$E, designSV, P = num_regress)
y_clean = cleaningY(v$E, designSV, P = num_regress)

## sanity checks on the regressed expression values
summary(y_clean['RBFOX3',]) # should be high
summary(y_clean['DRD1',]) # should be medium

###################################
# 2) grab the DEGs split by the DEGs to plot in neurons or glia
res = here(DATADIR, 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% readRDS()
path = here(PLOTDIR, 'tables','s3.2_table_numDEG_overlap_bycelltype.xlsx') 
res_overlap_list = excel_sheets(path) %>% set_names() %>%  map(read_excel, path = path)
sapply(res_overlap_list, nrow)

#############################################################
# 3) set up the heatmap plot parameters ahead of time
alpha = 0.05
gp_title = gpar(fontsize = 6, fontface = 'bold'); gp_label = gpar(fontsize = 3)
u0 = unit(0, "mm"); u1 = unit(1, "mm"); u3 = unit(2.4, "mm")
pal <- wes_palette("Zissou1", 21, type = "continuous")
pal = rainbow(12); pal = cividis(21)

row_params = list(title_gp = gpar(fontsize = 4, fontface = 'bold'), labels_gp = gp_label, 
                  column_gap = u1, grid_height = u1, grid_width = u3)
column_params = list(title_gp = gp_title, labels_gp = gp_label, ncol = 1, 
                     gap = unit(0, "cm"), row_gap = u0, column_gap = u0, 
                     grid_height = u1, grid_width = u1)

## reorder the samples to show cell type
df2 = df %>% arrange(celltype3, DSM.IV.OUD) %>% filter(celltype_class == 'Neuron')

## set up the column annotation plot parameters
column_ha1 = HeatmapAnnotation(
  'Cell type' = df2$celltype3, 'OUD Dx' = df2$DSM.IV.OUD,
  col = list(`Cell type` =  typecolors[4:9], `OUD Dx` = dx_col),
  simple_anno_size = u3, annotation_name_gp= gp_title,
  annotation_legend_param = column_params, show_legend = c(F, F))

to_plot= res$Neuron %>% filter(adj.P.Val.Between < alpha) %>% 
  filter(abs(logFC) > .2) %>% pull(gene) %>% unique()
mat1 = y_clean[to_plot,rownames(df2)]

## make the plot
plot_fn = here(PLOTDIR,'plots', 'figure3_celltypes_DEG_heatmap.neuron.pdf')
pdf(plot_fn, height = 80/in2mm, width = 130/in2mm, onefile = T)
ht1 = Heatmap(mat1, col = pal, name = "log2CPM", cluster_columns = FALSE, 
        show_column_names = F, row_names_gp = gp_label, 
        top_annotation = column_ha1, heatmap_legend_param = row_params)

ht1
dev.off()










to_plot2= res[names(typecolors)[4:9]] %>% rbindlist() %>% 
  filter(adj.P.Val.Between < 0.01) %>% group_by(gene) %>% 
  mutate(numCell = n(), varLfc = var(logFC, na.rm = T)) %>% 
  ungroup() %>% filter(numCell ==1, abs(logFC) > 1) %>%
  pull(gene) %>% unique()
to_plot2 = to_plot2[!to_plot2 %in% to_plot]
mat2 = y_clean[to_plot2,rownames(df2)]

ht2 = Heatmap(mat2, col = pal, name = "Subtype", cluster_columns = FALSE, 
              show_column_names = F, row_names_gp = gp_label,
              heatmap_legend_param = row_params, show_heatmap_legend = FALSE)

########################################################
## 4) plot all the glial DEGs w/ cell type overlaps
df2 = df %>% arrange(celltype3, DSM.IV.OUD) %>% filter(celltype_class == 'Glia')

## set up the column and row annotation plot parameters
column_ha2 = HeatmapAnnotation(
  'Cell type' = df2$celltype3, 'OUD Dx' = df2$DSM.IV.OUD,
  col = list(`Cell type` =  typecolors[10:15], `OUD Dx` = dx_col),
  simple_anno_size = u3, annotation_name_gp= gp_title,
  annotation_legend_param = column_params)

to_plot3= res$Glia %>% filter(adj.P.Val.Between < alpha) %>% 
  filter(abs(logFC) > 0.2) %>% pull(gene) %>% unique()
mat3 = y_clean[to_plot3,rownames(df2)]

to_plot4= res[names(typecolors)[10:15]] %>% rbindlist() %>% 
  filter(adj.P.Val.Between < 0.01) %>% group_by(gene) %>% 
  mutate(numCell = n(), varLfc = var(logFC, na.rm = T)) %>% ungroup() %>% 
  filter(numCell == 1, abs(logFC) > 0.2) %>% pull(gene) %>% unique()
to_plot4 = to_plot4[!to_plot4 %in% to_plot3]
mat4 = y_clean[to_plot4,rownames(df2)]

## make the plot
plot_fn = here(PLOTDIR,'plots', 's3.7_celltypes_DEG_heatmap.Glia.pdf')
pdf(plot_fn, height = 200/in2mm, width = 180/in2mm, onefile = T)
ht3 = Heatmap(mat3, col = pal,name = "log2CPM", cluster_columns = FALSE, 
              show_column_names = F, row_names_gp = gp_label, 
              top_annotation = column_ha2, heatmap_legend_param = row_params)

ht4 = Heatmap(mat4, col = pal, name = "Subtype", cluster_columns = FALSE, 
              show_column_names = F, row_names_gp = gp_label,
              heatmap_legend_param = row_params, show_heatmap_legend = FALSE)

ht3 %v% ht4
dev.off()




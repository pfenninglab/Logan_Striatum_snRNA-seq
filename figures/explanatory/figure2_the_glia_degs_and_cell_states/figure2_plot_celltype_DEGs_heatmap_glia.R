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
my_theme = theme_classic(base_size = 6)

###########################################
# 1) load the pseudobulk, regressed matrix
z_clean = readRDS(here(DATADIR, 'rdas', 'voomLimma_zscore_object_N222pb_regressedCovariate.rds'))

## grab the phenotype data
df = here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds') %>% 
  readRDS() %>% colData() %>% as.data.frame() %>% 
  dplyr::select(ID:CaseIdx) %>% 
  mutate(celltype3 = celltype3 %>% make.names %>% factor(names(typecolors)),
         celltype3 = droplevels(celltype3), 
         celltype_class = case_when(grepl('^D|^Int', celltype3)~ 'Neuron',  
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia')))
df = df[colnames(z_clean), ]

##################################################################
# 2) grab the DEGs split by the DEGs to plot in neurons or glia
res = here(DATADIR, 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% readRDS()
path = here('figures/explanatory/figure3_the_neuron_degs_and_cell_states',
            'tables','s3.2_table_numDEG_overlap_bycelltype.xlsx') 
res_overlap_list = path %>% excel_sheets() %>% set_names() %>% 
  map(read_excel, path = path) %>% lapply(function(x) x %>% filter(numCelltype >=2) )
sapply(res_overlap_list, nrow)

#############################################################
# 3) set up the heatmap plot parameters ahead of time
alpha = 0.05
gp_title = gpar(fontsize = 6, fontface = 'bold'); gp_label = gpar(fontsize = 4)
u0 = unit(0, "mm"); u2 = unit(2, "mm"); u3 = unit(3, "mm")
pal <- RColorBrewer::brewer.pal(11, 'PiYG')

row_params = list(title_gp = gp_title, labels_gp = gp_label, 
                  column_gap = u2, grid_height = u2, grid_width = u3, 
                  legend_height = u3)
column_params = list(title_gp = gp_title, labels_gp = gp_label, ncol = 1, 
                     gap = unit(0, "cm"), row_gap = u0, column_gap = u0, 
                     grid_height = u2, grid_width = u2)

## make the legends
lgd_zsc = Legend(col_fun = circlize::colorRamp2(seq(-4, 4,8/10 ), pal), 
                 title = "zscore", labels_gp = gp_label, gap = unit(0, "cm"), 
                 row_gap = u0, column_gap = u0, title_gp = gp_title,
                 grid_height = u2, grid_width = u2, legend_height = unit(6, 'mm'))
lgd_sex = Legend(labels = c('F', 'M'), title = "Sex", labels_gp = gp_label, 
                 legend_gp = gpar(fill =c('F' = 'black', 'M' = 'yellow')),
                 gap = unit(0, "cm"), row_gap = u0, column_gap = u0, legend_height = unit(6, 'mm'),
                 grid_height = u2, grid_width = u2, title_gp = gp_title)
lgd_oud = Legend(labels = names(dx_col), title = "OUD Dx", legend_height = unit(6, 'mm'),
                 legend_gp = gpar(fill = dx_col), labels_gp = gp_label,
                 gap = unit(0, "cm"), row_gap = u0, column_gap = u0, 
                 grid_height = u2, grid_width = u2, title_gp = gp_title)


########################################################
## 4) plot all the glial DEGs w/ cell type overlaps
df2 = df %>% arrange(celltype3, DSM.IV.OUD, Sex) %>% filter(celltype_class == 'Glia')%>% 
  dplyr::select(celltype3, DSM.IV.OUD, Sex)

## cell type legend
lgd_cell = Legend(labels = names(typecolors)[c(10:12,14:15)], title = "Cell type", 
                  legend_gp = gpar(fill =  typecolors[c(10:12,14:15)]), labels_gp = gp_label,
                  gap = unit(0, "cm"), row_gap = u0, column_gap = u0, legend_height = unit(6, 'mm'),
                  grid_height = u2, grid_width = u2, title_gp = gp_title)
lgd_list = packLegend(lgd_zsc, lgd_oud,  lgd_cell, lgd_sex, max_height = unit(80, 'mm'))

## set up the column and row annotation plot parameters
column_ha = HeatmapAnnotation(
  'Cell type' = df2$celltype3, 'OUD Dx' = df2$DSM.IV.OUD, 'Sex' = df2$Sex,
  col = list(`Cell type` =  typecolors[c(10:12,14:15)], `OUD Dx` = dx_col, 
             `Sex` = c('F' = 'black', 'M' = 'yellow')),
  simple_anno_size = u2, annotation_name_gp= gp_title,
  annotation_legend_param = column_params, show_legend = c(F, F, F))

## grab the genes to be plot on the heatmap
to_plot3= res$Glia %>% filter(adj.P.Val.Between < alpha) %>% 
  filter(abs(logFC) > 1) %>% pull(gene) %>% unique()
mat3 = z_clean[to_plot3,rownames(df2)]

bind_cols(df2, t(mat3)) %>% writexl::write_xlsx(here(PLOTDIR,'tables', 'figure2_celltypes_DEG_heatmap.glia.sourceData.xlsx'))

## grab the genes enriched in the pathways
pathways_genes = here(PLOTDIR,'tables', 'figure2_clustered_glia_gsea_pathway_network.xlsx') %>% 
  readxl::read_xlsx() %>% pull(leadingEdge) %>% sapply(str_split, ',') %>% unlist() %>% unique()

ind_mark = which(to_plot3 %in% pathways_genes)
row_ha = rowAnnotation(foo = anno_mark(at = ind_mark, labels = to_plot3[ind_mark], 
                                   link_width = unit(3, "mm"), labels_gp = gp_label))

## make the plot
plot_fn = here(PLOTDIR,'plots', 'figure2_celltypes_DEG_heatmap.Glia.pdf')
pdf(plot_fn, height = 70/in2mm, width = 121/in2mm, onefile = T)
ht2 = Heatmap(mat3, col = pal,name = "zscore", cluster_columns = T, 
        column_dend_height = u3, row_dend_width = u3,
        show_column_names = F, show_row_names = F,
        show_heatmap_legend = F, right_annotation = row_ha,
        top_annotation = column_ha, heatmap_legend_param = row_params)

draw(ht2, ht_gap = u0, annotation_legend_list = lgd_list)

dev.off()




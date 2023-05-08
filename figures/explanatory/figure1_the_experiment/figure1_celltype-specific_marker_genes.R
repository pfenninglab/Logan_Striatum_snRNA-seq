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

## down-sample oligos to same number of cells as astrocytes
ind_drop = which(obj_merged$celltype4 == 'Oligos') %>% sample(size = 54500)
obj_subset = obj_merged[, -ind_drop]
table(obj_subset$Case)
table(obj_subset$celltype4) %>% sort()

obj_subset = ScaleData(obj_subset)
Idents(obj_subset) = 'celltype4'
levels(obj_subset) = names(typecolors)[names(typecolors) %in% levels(obj_subset)]

######################################################################
# 2) compute the cell type marker genes w/ pseudobulking and DESeq2
diff_pb_gene_dir = here(PLOTDIR, 'plots', 'FindMarkersBulk_outs')
FindMarkersBulk( obj_subset, clus_ident = 'celltype4', sample_ident = 'Case',
  expfilt_counts = 1, expfilt_freq = 0.5, n_top_genes = 1000, pct.in = 25, alpha = 0.05,
  out_dir = diff_pb_gene_dir)

## find and read in the DESeq2 tables of one vs. all cell type marker genes
diff_fn = list.files(diff_pb_gene_dir, pattern = '_results.csv', full.names = T)
names(diff_fn) = basename(diff_fn) %>% ss('cluster_|_results', 2)
markers_df = diff_fn %>% lapply(read_csv, show_col_types = FALSE) %>% 
  bind_rows(.id = 'celltype')  %>% 
  dplyr::rename('gene' = '...1') %>% arrange(pvalue) %>% 
  ## readjust FDR across all cell type comparisons
  mutate(padj = lm_qvalue(pvalue, X=baseMean)$q, 
         pbonf = p.adjust(pvalue, 'bonferroni'),
         sig = case_when(padj < alpha ~ 'Significant', 
                         padj > alpha ~ 'Not significant')) %>% 
  relocate('pbonf', .after = 'padj')

markers_list = split(markers_df, markers_df$celltype)

## get the list of significant up-regulated marker genes by cell type
markers_signif = markers_list %>% lapply(function(df){
  df %>% filter(pbonf < alpha, log2FoldChange > 0) 
})

## filter down to marker genes that are not shared across cell types
markers_unique = markers_signif %>% lapply('[[', 'gene') %>% unlist() %>% table()
markers_unique = names(markers_unique)[markers_unique==1]
markers_signif = markers_signif %>% 
  lapply(function(df) df %>% filter(gene %in% markers_unique)) %>% 
  bind_rows() %>% arrange(pbonf)

markers_signif %>% dplyr::count(celltype) %>% arrange(desc(n))
markers_signif %>% dplyr::count(celltype) %>% pull(n) %>% summary()
markers_signif %>% dplyr::count(celltype) %>% pull(n) %>% sd()

markers_signif %>% filter(grepl('^D|Int', celltype)) %>% 
  dplyr::count(celltype) %>% pull(n) %>% summary()
markers_signif %>% filter(grepl('^D|Int', celltype)) %>% 
  dplyr::count(celltype) %>% pull(n) %>% sd()/sqrt(5)

markers_signif %>% filter(!grepl('^D|Int', celltype)) %>% 
  dplyr::count(celltype) %>% pull(n) %>% summary()
markers_signif %>% filter(!grepl('^D|Int', celltype)) %>% 
  dplyr::count(celltype) %>% pull(n) %>% sd()/sqrt(7)

##############################################################################
# 3) export the marker genes and all the markers computed with one vs. all pb
diff_markers_table = here(PLOTDIR, 'tables', 
                          'figure1_striatum_celltype_pseudobulk_marker_genes.xlsx')
markers_list2 = c(list('unique_signif_marker_genes' = markers_signif), markers_list) %>% 
  lapply(filter, padj < alpha)
markers_list2 %>% writexl::write_xlsx(diff_markers_table)

# # To load all sheets in a workbook, use lapply()
# markers_list2 = lapply(readxl::excel_sheets(diff_markers_table) %>% set_names(),
#                        readxl::read_excel, path = diff_markers_table)

#####################################################
# 4) plot the top marker genes across the cell types
markers_genes = markers_signif %>% group_by(celltype) %>% 
  filter(!grepl('\\.|^LINC|^TRAJ', gene)) %>% ## ignore genes w/ accession numbers
  top_n(5, log2FoldChange) %>% ungroup() %>% 
  mutate(celltype = celltype %>% factor(names(typecolors))) %>% 
  arrange(celltype) %>%  pull(gene)

# Aggregate across cluster-sample groups
groups <- obj_subset[[]][, c('celltype4', 'ID')] 
mat= GetAssayData(obj_subset, slot = 'scale.data', assay = 'RNA')[markers_genes,]
pb <- Matrix.utils::aggregate.Matrix(t(mat), groupings = groups, fun = "mean")
df = as.data.frame(pb) %>% rownames_to_column('groupings') %>% 
  mutate(celltype = gsub('_C-.*|_P-.*', '', groupings) %>% factor(names(typecolors))) %>% 
  pivot_longer(cols = -c(groupings, celltype), names_to = 'gene', values_to = 'value') %>% 
  arrange(celltype) %>% 
  mutate(gene = factor(gene, markers_genes), groupings = factor(groupings, unique(groupings)))

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



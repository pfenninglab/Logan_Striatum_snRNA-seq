ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(broom)
library(rcartocolor)
library(ggsci)

## main Seurat package snRNA-seq pacakges
library(SingleCellExperiment)
library(future)
library(here)

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
               'Oligos_Pre', 'Interneuron')
othertypes_col = c(carto_pal(4, "Safe"), 
                   mypal = pal_npg('nrc')(3),
                   carto_pal(length(othertypes) -7 , "Vivid"))
names(othertypes_col) = othertypes
othertypes_col = othertypes_col[-c(1:4)]
typecolors = c(subtypes_col, othertypes_col[c(10, 4:9)])

##############################################
# 1) load in the pseudo bulk TF regulon output 
save_pb_fn = here('figures/exploratory/AUCell_gene_set_activities', 'rdas',
                  'OUD_Striatum_pyscenicl_TF_modules.pseudoBulkByCelltype.rds')
aucell_pb = readRDS(save_pb_fn)
aucell_pb_list = split(aucell_pb, aucell_pb$TF)
table(aucell_pb$celltype3)
typecolors = typecolors[names(typecolors) %in% aucell_pb$celltype3]

############################################################################
# 2) compute the cell type positive TF regulator w/ pseudobulking and DESeq2
aucell_diff_df = lapply(set_names(names(typecolors)), function(ct){
  ## loop through each TF to calculate differential TF-modules per cell type
  diff_df = parallel::mclapply(aucell_pb_list, function(df){
    if(sum(df$celltype3 == ct) >4 ){
      df = df %>% mutate(
        celltype = ifelse(celltype3 == ct, 'Target', 'Not Target')
      )
      if(length(unique(df$celltype)) ==2 ){
      out = lm(AUCell ~ celltype + Age + Sex + PMI + Region, data = df) %>% 
        tidy() %>%  filter(term == 'celltypeTarget')
      return(out)
    }}
  }, mc.cores = 16) %>% bind_rows(.id = 'TF_module')
  return(diff_df)
}) %>% bind_rows(.id = 'celltype')

## find and read in the DESeq2 tables of one vs. all cell type marker genes
markers_df = aucell_diff_df  %>% arrange(p.value) %>% dplyr::select(-term) %>% 
  ## readjust FDR across all cell type comparisons
  mutate(padj = p.adjust(p.value, 'fdr'), 
         pbonf = p.adjust(p.value, 'bonferroni'),
         sig = case_when(padj < alpha ~ 'Significant', 
                         padj > alpha ~ 'Not significant')) %>% 
  relocate('pbonf', .after = 'padj')

markers_list = split(markers_df, markers_df$celltype)

## get the list of significant up-regulated marker genes by cell type
markers_signif = markers_list %>% lapply(function(df){
  df %>% filter(padj < alpha, estimate > 0) 
})

markers_signif %>% bind_rows() %>% dplyr::count(celltype)

## filter down to marker genes that are not shared across cell types
markers_unique = markers_signif %>% lapply('[[', 'TF_module') %>% unlist() %>% table()
markers_unique = names(markers_unique)[markers_unique<=3]
markers_signif = markers_signif %>% 
  lapply(function(df) df %>% filter(TF_module %in% markers_unique)) %>% 
  bind_rows() %>% arrange(padj)

markers_signif %>% dplyr::count(celltype)
markers_signif %>% dplyr::count(celltype) %>% arrange(desc(n))

markers_signif %>% dplyr::count(celltype) %>% pull(n) %>% summary()
markers_signif %>% dplyr::count(celltype) %>% pull(n) %>% sd()/sqrt(12)

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
                          'figure1_striatum_celltype_pseudobulk_marker_TF-modules.xlsx')
markers_list2 = c(list('unique_signif_marker_TF-modules' = markers_signif), markers_list) %>% 
  lapply(filter, padj < alpha)
markers_list2 %>% writexl::write_xlsx(diff_markers_table)

# # To load all sheets in a workbook, use lapply()
# markers_list2 = lapply(readxl::excel_sheets(diff_markers_table) %>% set_names(),
#                        readxl::read_excel, path = diff_markers_table)


#####################################################
# 4) plot the top marker genes across the cell types
markers_TF_modules = markers_signif %>% group_by(celltype) %>% 
  top_n(7, estimate) %>% ungroup() %>% 
  mutate(celltype = celltype %>% factor(names(typecolors))) %>% 
  arrange(celltype) %>% pull(TF_module) %>% unique()

# Aggregate across cluster-sample groups
df = aucell_pb %>% filter(TF %in% markers_TF_modules) %>% 
  group_by(TF) %>% 
  mutate(AUCellz = (AUCell - min(AUCell, na.rm = T))/sd(AUCell)) %>% 
  mutate(celltype = factor(celltype3, names(typecolors))) %>% 
  arrange(celltype) %>% 
  mutate(TF_module = factor(TF, markers_TF_modules), 
         groupings = paste(celltype, ID))

df2 = df %>% dplyr::select(c(groupings, TF_module, celltype, AUCellz)) %>% 
  pivot_wider(id_cols = c(groupings, celltype), names_from = TF_module, 
              values_from = AUCellz, values_fill = 0) %>% 
  pivot_longer(-c(groupings, celltype), values_to = 'AUCellz', names_to = 'TF_module') %>% 
  mutate(TF_module = factor(TF_module, markers_TF_modules))

## make marker TF_module plot
plot_fn = here(PLOTDIR, 'plots', 'figure1_pseudobulk_celltype_markerTF_module_heatmap.pdf')
pdf(plot_fn, height = 60/in2mm, width = 90/in2mm, onefile = T)
ggplot(df2, aes(x = groupings, y = TF_module, fill = AUCellz))+
  geom_tile() + 
  scale_fill_distiller(palette = "Spectral", 'Z-score') + 
  scale_y_discrete(limits=rev) + 
  facet_grid(~celltype, scales = 'free_x') +
  my_theme + xlab('Biological Sample') + 
  ylab('Transcription factor-gene regulatory module') + 
  theme_classic(base_size = 4) +
  theme( legend.key.height = unit(2, 'mm'), legend.position = 'bottom', 
         axis.text.x = element_blank(), axis.ticks
         = element_blank(), 
         strip.text.x = element_text(size = 3),
         legend.box.margin=margin(-2,-2,-2,-2))
dev.off()



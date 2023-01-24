ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(data.table)
library(ggsci)
library(here)
library(readxl)
library(Matrix.utils)

DATADIR='data/tidy_data/differential_expression_analysis'

## make for this subdirs
PLOTDIR='figures/explanatory/figure5_sex_interaction_degs'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

###################################
# 0) pre-set plotting aesthetics
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes %>% make.names()

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 'All', 'Neuron', 'Glia',
               'Astrocytes', 'Endothelial', 'Microglia', 'Mural', 'Oligos', 
               'Oligos_Pre', 'Interneuron')
othertypes_col = c(carto_pal(4, "Safe"), 
                   mypal = pal_npg('nrc')(3),
                   carto_pal(length(othertypes) -7 , "Vivid"))
names(othertypes_col) = othertypes
othertypes_col = othertypes_col[-c(1:4)]
typecolors = c(othertypes_col[1:3], subtypes_col, othertypes_col[c(10, 4:9)])

in2mm<-25.4
my_theme = theme_classic(base_size = 5)

#########################################################
# 1) grab the phenotype data & the cleaned gene log2 matrix
y_clean = readRDS(here(DATADIR, 'rdas', 'voomLimma_norm_object_N222pb_regressedCovariate.rds'))

df = here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds') %>% 
  readRDS() %>% SummarizedExperiment::colData() %>% as.data.frame() %>% 
  dplyr::select(ID:CaseIdx) %>% 
  mutate(celltype3 = celltype3 %>% make.names %>% factor(names(typecolors)),
         celltype3 = droplevels(celltype3), 
         celltype_class = case_when(grepl('^D|^Int', celltype3)~ 'Neuron',  
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia')))

## reorder
df = df[colnames(y_clean), ] 

## sanity checks on the regressed expression values
summary(y_clean['RBFOX3',]) # should be high
summary(y_clean['DRD1',]) # should be medium
all(colnames(y_clean) == rownames(df)) # should be true

###################################
# 2) grab the DEGs split by the DEGs to plot in neurons or glia
alpha = 0.05

path = here(PLOTDIR, 'tables','s5.2_table_numDEG_overlap_bySexInteraction.xlsx')
res_overlap_list = path %>% excel_sheets() %>% set_names() %>% 
  map(read_excel, path = path) %>% lapply(function(x) x %>% filter(numCelltype >=3) )
sapply(res_overlap_list, nrow)

res = here(DATADIR, 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds') %>% readRDS()

###################################
# 3) compute average log2 CPM for plotting the degree significance
y_summary = data.frame(celltype = df$celltype3, t(y_clean) %>% as.data.frame()) %>% 
  pivot_longer(cols = -celltype, names_to = 'gene', values_to = 'Y') %>% 
  group_by(celltype, gene) %>% summarise(Y_max = mean(Y), Y_sd = sd(Y)) %>% 
  mutate(match = paste0(celltype,'.',gene)) %>% column_to_rownames('match')

res_signif = res %>% rbindlist() %>% dplyr::select(-contains('dir')) %>% 
  pivot_longer(cols = -c(celltype, gene), names_sep = '_Sex', names_to = c('name', 'Sex')) %>% 
  pivot_wider(id_cols =c(celltype, gene, Sex)) %>% 
  filter(adj.P.Val.Between < alpha) %>% 
  filter(!celltype %in% c('All', 'Neuron', 'Glia')) %>% 
  dplyr::rename( 'celltype3' = 'celltype') %>% 
  mutate(annotation = case_when(adj.P.Val.Between < 0.001 ~ '***', 
                                adj.P.Val.Between < 0.01 ~ '**', 
                                adj.P.Val.Between < alpha ~ '*'), 
         Y_max = y_summary[paste0(celltype3,'.',gene), 'Y_max']  + 
           2 * (y_summary[paste0(celltype3,'.',gene), 'Y_sd']),
         celltype_class = case_when(grepl('^D|^Int', celltype3)~ 'Neuron',  
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia'))) %>% 
  filter(!is.na(Y_max))

#############################################################
## 4) make the plots of the log2CPMs by cell type and OUD dx
dir.create(here(PLOTDIR,'plots','deg_plots'), showWarnings = F)

# plot all the neuronal DEGs w/ >=2 cell type overlaps
to_plot= c('FKBP5', 'SLC17A6')
for(g in to_plot){
  plot_fn = here(PLOTDIR,'plots', paste0('figure5.sexInteraction_celltypes_DEG_boxPlot.',g,'.pdf'))
  pdf(plot_fn, height = 40/in2mm, width = 55/in2mm)
  df2 = cbind(df, Y = y_clean[g,rownames(df)])
  p = ggplot(df2, aes(x = celltype3, y = Y)) + 
    geom_boxplot(outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
    geom_point( alpha = .7, size = .5, position = position_jitterdodge(),
                aes(color = DSM.IV.OUD)) + 
    geom_text(data =  res_signif %>% filter(gene %in% g),
              size = 2, hjust = 'inside', aes(y = Y_max, label = annotation, Sex = Sex)) +
    facet_grid(Sex ~ celltype_class, scales = 'free', space = 'free') + 
    my_theme + scale_fill_npg() + 
    ylab(paste(g, 'log2(CPM) | SVs')) + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
    theme(legend.position = 'none', 
          legend.spacing.y = unit(2, 'cm'), 
          legend.key.size = unit(.5, "cm")) 
  print(p)
  dev.off()
}



#############################################################
## 5) make the supplemental plots of the log2CPMs by cell type and OUD dx


# plot all the neuronal DEGs w/ >=2 cell type overlaps
to_plot= res_overlap_list[['All_overlap']] %>% pull(gene)
plot_fn = here(PLOTDIR,'plots', 's5.4_celltypes_DEG_boxPlot.All.pdf')
pdf(plot_fn, height = 250/in2mm, width = 180/in2mm, onefile = T)
for(g in split(sort(to_plot), ceiling(seq_along(to_plot)/8))){
  df2 = df %>% cbind( t(y_clean[g,rownames(df)]) ) %>% 
    pivot_longer(cols = all_of(g), names_to = 'gene', values_to = 'Y')
  p = ggplot(df2, aes(x = celltype3, y = Y)) + 
    geom_boxplot(outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
    geom_point( alpha = .7, size = 1, position = position_jitterdodge(),
                aes(color = DSM.IV.OUD, shape = Sex)) + 
    geom_text(data =  res_signif %>% filter(gene %in% g),
              size = 2, hjust = 'inside', aes(y = Y_max, label = annotation, Sex = Sex)) +
    facet_wrap(~gene + Sex, scales = 'free', ncol = 8) + my_theme +  scale_fill_npg() + 
    ylab(expression(log[2]*'(CPM)|Covariates & SVs')) + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
    theme(legend.position = 'bottom', 
          legend.spacing.y = unit(2, 'cm'), 
          legend.key.size = unit(.5, "cm")) 
  print(p)
}
dev.off()


# plot all the neuronal DEGs w/ >=2 cell type overlaps
to_plot= res_overlap_list[['Neuron_overlap']] %>% pull(gene)
plot_fn = here(PLOTDIR,'plots', 's5.5_celltypes_DEG_boxPlot.Neuron.pdf')
pdf(plot_fn, height = 250/in2mm, width = 180/in2mm, onefile = T)
for(g in split(sort(to_plot), ceiling(seq_along(to_plot)/12))){
  df2 = df %>% cbind( t(y_clean[g,rownames(df)]) ) %>% 
    pivot_longer(cols = all_of(g), names_to = 'gene', values_to = 'Y') %>% 
    filter(celltype_class == 'Neuron')
  p = ggplot(df2, aes(x = celltype3, y = Y)) + 
    geom_boxplot(outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
    geom_point( alpha = .7, size = 1, position = position_jitterdodge(),
                aes(color = DSM.IV.OUD, shape = Sex)) + 
    geom_text(data =  res_signif %>% filter(gene %in% g) %>% 
                filter(celltype_class == 'Neuron'), size = 2, hjust = 'inside',
              aes(y = Y_max, label = annotation)) +
    facet_wrap(~gene + Sex, scales = 'free', ncol = 4) + 
    my_theme +  scale_fill_npg() + 
    ylab(expression(log[2]*'(CPM)|Covariates & SVs')) + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
    theme(legend.position = 'bottom', 
          legend.spacing.y = unit(2, 'cm'), 
          legend.key.size = unit(.5, "cm")) 
  print(p)
}
dev.off()


# plot all the neuronal DEGs w/ >=2 cell type overlaps
to_plot= res_overlap_list[['Glia_overlap']] %>% pull(gene)
plot_fn = here(PLOTDIR,'plots', 's5.6_celltypes_DEG_boxPlot.Glia.pdf')
pdf(plot_fn, height = 200/in2mm, width = 180/in2mm, onefile = T)
for(g in split(sort(to_plot), ceiling(seq_along(to_plot)/12))){
  df2 = df %>% cbind( t(y_clean[g,rownames(df)]) ) %>% 
    pivot_longer(cols = all_of(g), names_to = 'gene', values_to = 'Y') %>% 
    filter(celltype_class == 'Glia')
  p = ggplot(df2, aes(x = celltype3, y = Y)) + 
    geom_boxplot(outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
    geom_point( alpha = .7, size = 1, position = position_jitterdodge(),
                aes(color = DSM.IV.OUD, shape = Sex)) + 
    geom_text(data =  res_signif %>% filter(gene %in% g) %>% 
                filter(celltype_class == 'Glia'), size = 2, hjust = 'inside',
              aes(y = Y_max, label = annotation)) +
    facet_wrap(~gene + Sex, scales = 'free', ncol = 4) + 
    my_theme +  scale_fill_npg() + 
    ylab(expression(log[2]*'(CPM)|Covariates & SVs')) + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))+
    theme(legend.position = 'bottom', 
          legend.spacing.y = unit(2, 'cm'), 
          legend.key.size = unit(.5, "cm")) 
  print(p)
}
dev.off()




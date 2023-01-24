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
library(SingleCellExperiment)
library(Matrix.utils)

DATADIR='data/tidy_data/differential_expression_analysis'

## make for this subdirs
PLOTDIR='figures/exploratory/compare_to_neural_activity_mouse'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)
in2mm<-25.4
my_theme = theme_classic(base_size = 5)

###########################################
# 1) load the pseudobulk, regressed matrix
y_clean = readRDS(here(DATADIR, 'rdas', 'voomLimma_norm_object_N222pb_regressedCovariate.rds'))

## grab the phenotype data
df = here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds') %>% 
  readRDS() %>% colData() %>% as.data.frame() %>% 
  dplyr::select(ID:CaseIdx) %>% 
  mutate(celltype3 = celltype3 %>% make.names %>% factor(names(typecolors)),
         celltype3 = droplevels(celltype3), 
         celltype_class = case_when(grepl('^D|^Int', celltype3)~ 'Neuron',  
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia')))
df = df[colnames(y_clean), ]

###################################
# 2) grab the DEGs split by the DEGs to plot in neurons or glia
alpha = 0.05
res = here(DATADIR, 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% readRDS()

######################################################################
# 3) compute average log2 CPM for plotting the degree significance
y_summary = data.frame(celltype = df$celltype3, t(y_clean) %>% as.data.frame()) %>% 
  pivot_longer(cols = -celltype, names_to = 'gene', values_to = 'Y') %>% 
  group_by(celltype, gene) %>% summarise(Y_max = mean(Y), Y_sd = sd(Y)) %>% 
  mutate(match = paste0(celltype,'.',gene)) %>% column_to_rownames('match')

res_signif = res %>% rbindlist() %>% filter(adj.P.Val.Between < alpha) %>% 
  filter(!celltype %in% c('All', 'Neuron', 'Glia')) %>% 
  dplyr::rename( 'celltype3' = 'celltype') %>% 
  mutate(annotation = case_when(adj.P.Val.Between < 0.001 ~ '***', 
                                adj.P.Val.Between < 0.01 ~ '**', 
                                adj.P.Val.Between < alpha ~ '*'), 
         Y_max = y_summary[paste0(celltype3,'.',gene), 'Y_max']  + 
           2.5 * (y_summary[paste0(celltype3,'.',gene), 'Y_sd']),
         celltype_class = case_when(grepl('^D|^Int', celltype3)~ 'Neuron',  
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia')))


#############################################################
## 4) make the plots of the log2CPMs by cell type and OUD dx
dir.create(here(PLOTDIR,'plots','deg_plots'), showWarnings = F)

to_plot= c('ARC','FOS','FOSB', 'EGR1','EGR2','EGR3', 'NR4A1')
to_plot = to_plot[to_plot%in% rownames(y_clean)]
plot_fn = here(PLOTDIR,'plots', 'compare_neural_activity_IEG_boxplot.pdf')
pdf(plot_fn, height = 70/in2mm, width = 60/in2mm, onefile = T)
df2 = df %>% cbind( t(y_clean[to_plot,rownames(df)]) ) %>% 
  pivot_longer(cols = all_of(to_plot), names_to = 'gene', values_to = 'Y')
p = ggplot(df2, aes(x = celltype3, y = Y)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
  geom_point( alpha = .7, size = .7, position = position_jitterdodge(),
              aes(color = DSM.IV.OUD, shape = Sex)) + 
  geom_text(data =  res_signif %>% filter(gene %in% to_plot),
            size = 2, vjust = 'inside', aes(y = Y_max, label = annotation)) +
  facet_grid(gene~celltype_class, scales = 'free', space = 'free_x') + 
  my_theme +  scale_fill_npg() + 
  ylab(expression(log[2]*'(CPM)|Covariates & SVs')) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))+
  theme(legend.position = 'bottom', 
        legend.spacing.y = unit(0, 'cm'), 
        legend.key.size = unit(1, "mm"), 
        legend.margin=margin(-1,0,0,0, "mm")) 
print(p)
dev.off()


to_plot= c('NPAS4', 'BDNF')
to_plot = to_plot[to_plot%in% rownames(y_clean)]
plot_fn = here(PLOTDIR,'plots', 'compare_neural_activity_SRG_boxplot.pdf')
pdf(plot_fn, height = 70/in2mm, width = 60/in2mm, onefile = T)
df2 = df %>% cbind( t(y_clean[to_plot,rownames(df)]) ) %>% 
  pivot_longer(cols = all_of(to_plot), names_to = 'gene', values_to = 'Y')
p = ggplot(df2, aes(x = celltype3, y = Y)) + 
  geom_boxplot(outlier.shape = NA, aes(fill = DSM.IV.OUD)) + 
  geom_point( alpha = .7, size = .7, position = position_jitterdodge(),
              aes(color = DSM.IV.OUD, shape = Sex)) + 
  geom_text(data =  res_signif %>% filter(gene %in% to_plot),
            size = 2, vjust = 'inside', aes(y = Y_max, label = annotation)) +
  facet_grid(gene~celltype_class, scales = 'free', space = 'free_x') + 
  my_theme +  scale_fill_npg() + 
  ylab(expression(log[2]*'(CPM)|Covariates & SVs')) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 20, vjust = 1, hjust=1))+
  theme(legend.position = 'bottom', 
        legend.spacing.y = unit(0, 'cm'), 
        legend.key.size = unit(1, "mm"), 
        legend.margin=margin(-1,0,0,0, "mm")) 
print(p)
dev.off()



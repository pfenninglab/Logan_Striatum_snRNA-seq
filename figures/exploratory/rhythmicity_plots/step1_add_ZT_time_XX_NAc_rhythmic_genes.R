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
library(tidyverse)
library(data.table)
library(tidymodels)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

DATADIR='data/tidy_data/differential_expression_analysis'
PLOTDIR = 'figures/exploratory/rhythmicity_plots'

sapply(here(PLOTDIR, c('plots', 'tables', 'rdas')), dir.create, 
       recursive = T, showWarnings = F)

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
z_clean = readRDS(here(DATADIR, 'rdas', 'voomLimma_norm_object_N222pb.rds'))

zt_case = readxl::read_xlsx(here(PLOTDIR, 'tables', 'Striatum_ZT.xlsx')) %>% 
  mutate(Case= as.character(Case)) %>% deframe()

## grab the phenotype data
df = here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds') %>% 
  readRDS() %>% colData() %>% as.data.frame() %>% 
  dplyr::select(ID:CaseIdx) %>% 
  mutate(ZT = zt_case[as.character(Case)],
    celltype3 = celltype3 %>% make.names %>% factor(names(typecolors)),
         celltype3 = droplevels(celltype3), 
         celltype_class = case_when(grepl('^D|^Int', celltype3)~ 'Neuron',  
                                    TRUE ~ 'Glia'), 
         celltype_class = factor(celltype_class , c('Neuron', 'Glia')))
df = df[colnames(z_clean), ]

table(df$Case %in% names(zt_case))
table(df$ZT)

rhythm_str = read_csv(here('data/tidy_data/Xue2022_Rhythmicity_tables/NAC_OUDvsCONT.csv')) %>% 
  filter(R2losePvalue == 'TRUE' | lossSig == 'TRUE') %>% 
  dplyr::select(-c(R2losePvalue:intshiftPvalue)) %>% 
  pivot_longer(cols = gainSig:lossSig) %>% 
  filter(value == 'TRUE', genes %in% rownames(z_clean)) %>% 
  dplyr::select(-c(`...1`, 'value')) %>% 
  split(x = .$genes, f = .$name)


##################################
# 2) make the plots for each gene
label = names(rhythm_str)
to_plot = rhythm_str[['gainSig']] 

for (label in  names(rhythm_str)) {
  to_plot = rhythm_str[[label]] 
  for(gene in to_plot){
    
  df2 =df %>% mutate(
    y = z_clean[gene,rownames(df)] %>% as.matrix() %>% as.numeric(),
    celltype = celltype3,
    Dx = case_when(DSM.IV.OUD== 'CTL' ~ 'UC', T ~ 'OUD'),
  ) %>% 
    nest(data = -c(DSM.IV.OUD, celltype3)) %>% 
    filter()
  
  df3 = lapply(df2$data, function(d){
    ## fit the sine curve to get amplitude and phase
    fit <- lm(y ~ sin(pi/12*ZT)+cos(pi/12*ZT), data = d)
    d2 = data.frame(ZT = seq(-1:17)) %>% 
      mutate(celltype =unique( d$celltype), 
             Dx = unique(d$Dx))
    ## get the sine curve fits
    d2$y_pred <- predict(fit, newdata = d2)
    return(d2)
  }) %>% rbindlist()
  
  df2 = df2 %>% unnest(data)
  
  out_fn = here(PLOTDIR, 'plots', paste0('snRNA_', label, '.', gene,'.pdf'))
  pdf(out_fn, height = 3, width = 6 )
  # Assuming you have a data frame called 'df' with x and y values
  p= ggplot(df2, aes(x = ZT, y = y, fill = Dx)) +
    geom_point(pch = 21) + 
    geom_line(data = df3, aes(x = ZT, y = y_pred, color = Dx))+ 
    facet_wrap(~celltype, scales = 'free_y', nrow = 2 ) + 
    theme_bw(base_size = 8) + 
    theme(legend.position = 'bottom') + 
    ylab(gene) + xlim(c(-3,18))
  print(p)
  dev.off()
}}



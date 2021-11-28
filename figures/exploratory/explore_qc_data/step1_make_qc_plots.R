#
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(tidyverse)
library(here)
library(rcartocolor)

DATADIR='data/raw_data'
PROJDIR='figures/exploratory/explore_qc_data'

## make data dirs
here(PROJDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

## read in the phenotype data
pd = read_csv(here(DATADIR, 'tables', 'OUD1_snRNA_seq_sampleSheet.csv')) %>%
  rename_with(make.names) %>%
  mutate(Center = ifelse(Region =='DLPFC', 'Singulomics', 'Boston University'), 
         match1 = Sample.ID, match2 = paste(Sample.ID, Region, sep = '.'),
         Sample.ID = gsub('HU', 'D-', Sample.ID), 
         Subject = paste0('HU', ss(Sample.ID, '-', 2)))

#############################################
## read in the STARsolo mapping QC data tables
map_qc_df = hciR::read_STAR(here(DATADIR, 'STARsolo_out')) %>%
  dplyr::rename( 'match1' =  'sample') %>%
  mutate(stat = factor(stat, unique(stat)))  

map_qc_df2 = inner_join(x = pd %>% dplyr::select(c(Sample.ID, match1, match2, Region, Pair, Case, Center)),
                        y = map_qc_df, by = 'match1') %>%
  dplyr::select(-c(match2, match1))

plot_fn = here(PROJDIR, 'plots', 'Logan-OUD-snRNAseq.mapping_QC_stats.BUrun1_Singulomics.pdf')
pdf(plot_fn, width = 6, height =4)
for(x in levels(map_qc_df$stat)){
  tmp = map_qc_df2 %>%filter(stat %in% x)
  if(! all(tmp$value ==0)){
  pp= ggplot(data = tmp, aes(x = Sample.ID, y = value, fill = Region)) +
    geom_bar(stat = 'identity') +
    scale_fill_carto_d(name = 'Tissue', palette = 'Safe') + ylab(x) + 
    facet_grid(~Center + Region, space = 'free', scales = 'free_x') + 
    theme_bw(base_size = 10) + 
    theme(legend.position = 'bottom')
  print(pp)
  }
}
dev.off()

tab_fn = here(PROJDIR, 'tables', 'Logan-OUD-snRNAseq.mapping_QC_stats.BUrun1_Singulomics.xlsx')
map_qc_df2 %>%pivot_wider(values_from = 'value', names_from = 'stat') %>%
  writexl::write_xlsx(tab_fn)






################################################################
######## read in STAR solo Count QC
summary_qc_fn = list.files(path = here(DATADIR, 'STARsolo_out'), recursive = T,
                       pattern = 'Summary.csv', full.names = T)
names(summary_qc_fn) = ss(summary_qc_fn, '/', 9) %>% ss('\\.Solo',1)
sumary_qc_df = summary_qc_fn %>% 
  lapply(read_csv, col_names = c('stat', 'value'),show_col_types = FALSE) %>%
  data.table::rbindlist(idcol = 'match2') %>%
  mutate(stat = factor(stat, unique(stat)))

sumary_qc_df2 = inner_join(x = pd%>% dplyr::select(c(Sample.ID, match1, match2, Region, Pair, Case, Center)),
                           y = sumary_qc_df, by = 'match2') %>%
  dplyr::select(-c(match2, match1))

plot_fn2 = here(PROJDIR, 'plots', 'Logan-OUD-snRNAseq.summary_QC_stats.BUrun1_Singulomics.pdf')
pdf(plot_fn2, width = 6, height =4)
for(x in levels(sumary_qc_df$stat)){
  tmp = sumary_qc_df2 %>%filter(stat %in% x)
  if(! all(tmp$value ==0)){
    pp= ggplot(data = tmp, aes(x = Sample.ID, y = value, fill = Region)) +
      geom_bar(stat = 'identity') +
      scale_fill_carto_d(name = 'Tissue', palette = 'Safe') + ylab(x) + 
      facet_grid(~Center + Region, space = 'free', scales = 'free_x') + 
      theme_bw(base_size = 10) + 
      theme(legend.position = 'bottom')
    print(pp)
  }
}
dev.off()

tab_fn2 = here(PROJDIR, 'tables', 'Logan-OUD-snRNAseq.summary_QC_stats.BUrun1_Singulomics.xlsx')
sumary_qc_df2 %>%pivot_wider(values_from = 'value', names_from = 'stat') %>%
  writexl::write_xlsx(tab_fn2)



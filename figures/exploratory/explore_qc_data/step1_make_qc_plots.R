#
ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

library(tidyverse)
library(here)
library(rcartocolor)
# devtools::install_github("HuntsmanCancerInstitute/hciR")

DATADIR='data/raw_data'
PROJDIR='figures/exploratory/explore_qc_data'

## make data dirs
here(PROJDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

## read in the phenotype data
pd = here('data/tidy_data', 'tables', 'LR_RM_BUMSR_Library_Index_Key_Run1_Run2.csv') %>%
  read_csv() %>%
  rename_with(make.names) %>%
  mutate(Region = case_when(grepl('P',Sample) ~ 'Putamen',
                            grepl('C',Sample) ~ 'Caudate'),
         Subject = paste0('HU', ss(Sample, '-', 2))) %>%
  dplyr::select(-c(Lane, Subject))


#############################################
## read in the STARsolo mapping QC data tables
map_qc_df = hciR::read_STAR(here(DATADIR, 'STARsolo_out')) %>%
  dplyr::rename( 'Sample' =  'sample') %>%
  mutate(stat = factor(stat, unique(stat)))  
table(map_qc_df$Sample[!map_qc_df$Sample %in% pd$Sample])

map_qc_df2 = inner_join(x = pd %>% dplyr::select(c(Sample, Region)),
                        y = map_qc_df) 

plot_fn = here(PROJDIR, 'plots', 'Logan-OUD-snRNAseq.mapping_QC_stats.BU_Run0-7.pdf')
pdf(plot_fn, width = 6, height =4)
for(x in levels(map_qc_df$stat)){
  tmp = map_qc_df2 %>%filter(stat %in% x)
  if(! all(tmp$value ==0)){
  pp= ggplot(data = tmp, aes(x = Sample, y = value, fill = Region)) +
    geom_bar(stat = 'identity') +
    scale_fill_carto_d(name = 'Tissue', palette = 'Safe') + ylab(x) + 
    facet_grid(~ Region, space = 'free', scales = 'free_x') + 
    theme_bw(base_size = 10) + 
    theme(legend.position = 'bottom') + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  print(pp)
  }
}
dev.off()


tab_fn = here(PROJDIR, 'tables', 'Logan-OUD-snRNAseq.mapping_QC_stats.BU_Run0-7.xlsx')
map_qc_df2 %>%pivot_wider(values_from = 'value', names_from = 'stat') %>%
  writexl::write_xlsx(tab_fn)






################################################################
######## read in STAR solo Count QC
summary_qc_fn = list.files(path = here(DATADIR, 'STARsolo_out'), recursive = T,
                       pattern = 'Summary.csv', full.names = T) %>%
  stringr::str_subset(pattern = 'GeneFull')
names(summary_qc_fn) = ss(summary_qc_fn, '/', 9) %>% ss('\\.Solo',1)
sumary_qc_df = summary_qc_fn %>% 
  lapply(read_csv, col_names = c('stat', 'value'),show_col_types = FALSE) %>%
  data.table::rbindlist(idcol = 'Sample') %>%
  mutate(stat = factor(stat, unique(stat)))

sumary_qc_df2 = inner_join(x = pd%>% dplyr::select(c(Sample, Region)),
                           y = sumary_qc_df, by = 'Sample') 


plot_fn2 = here(PROJDIR, 'plots', 'Logan-OUD-snRNAseq.summary_QC_stats.BU_Run0-7.pdf')
pdf(plot_fn2, width = 6, height =4)
for(x in levels(sumary_qc_df$stat)){
  tmp = sumary_qc_df2 %>%filter(stat %in% x)
  if(! all(tmp$value ==0)){
    pp= ggplot(data = tmp, aes(x = Sample, y = value, fill = Region)) +
      geom_bar(stat = 'identity') +
      scale_fill_carto_d(name = 'Tissue', palette = 'Safe') + ylab(x) + 
      facet_grid(~ Region, space = 'free', scales = 'free_x') + 
      theme_bw(base_size = 10) + 
      theme(legend.position = 'bottom') + 
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    print(pp)
  }
}
dev.off()

tab_fn2 = here(PROJDIR, 'tables', 'Logan-OUD-snRNAseq.summary_QC_stats.BU_Run0-7.xlsx')
sumary_qc_df2 %>%pivot_wider(values_from = 'value', names_from = 'stat') %>%
  writexl::write_xlsx(tab_fn2)



library(tidyverse)
library(data.table)
library(rtracklayer)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/ldsc_interferon_response'
dir.create(here(DATADIR, 'bed'), showWarnings = F, recursive = T)
dir.create(here(DATADIR, 'rdas'), showWarnings = F, recursive = T)
dir.create(here(DATADIR, 'tables'), showWarnings = F, recursive = T)

########################################
### 1) import the list of peaks by group
df = readxl::read_xlsx(here(DATADIR, 'monocyte_interferon_peaks.xlsx')) %>% 
  as.data.table() %>% rename_with(make.names)

df_bg = df %>% filter(Group == 'bg_cell_types')

## pivot to group cell types by their matched Naive files
df_ifn = df %>% filter(Group == 'interferon_test') %>% 
  pivot_wider(id_cols = c(Group, Source, Celltype), names_from = 'Treatment', 
               values_from = 'File.Name') %>% 
  pivot_longer(cols = c(IFNb, IFNg, `IFNg-LPS`), names_to = 'Treatment', 
               values_to = 'File.Tx') %>% 
  dplyr::rename('File.Naive' = 'Naive') %>% filter(lengths(File.Tx)>=1) %>% 
  mutate(Source = str_replace(Source, '\\.,', ''), 
         Label.Tx = paste(Source, Celltype, Treatment, 'Unique', sep = '_'), 
         Label.Tx = make.names(Label.Tx), 
         Label.Shared = paste(Source, Celltype, Treatment, 'SharedWithNaive', sep = '_'), 
         Label.Shared = make.names(Label.Shared))

## save the tables
df_ifn %>% saveRDS(here(DATADIR, 'rdas', 'interferon_macrophages_monocyte_annotations.rds'))
df_ifn %>% write_csv(here(DATADIR, 'tables', 'interferon_macrophages_monocyte_annotations.csv'))
df %>% saveRDS(here(DATADIR, 'rdas', 'blood_celltypes_and_interferon_annotations.rds'))

####################################
## 2) read in the peaks as GRanges
peaks_fn = list.files(here(DATADIR,'All'), full.names = T)
names(peaks_fn) = basename(peaks_fn) %>% ss('\\.')
peaks_fn = peaks_fn[df$File.Name]
peaks = lapply(peaks_fn, import)

## gather and reduce to one set of treatment and naive cell type peaks
naive_peaks = lapply(seq(nrow(df_ifn)), function(i){
  peaks[df_ifn$File.Naive[[i]]] %>% GRangesList() %>% unlist() %>% reduce()
})

treat_peaks = lapply(seq(nrow(df_ifn)), function(i){
  peaks[df_ifn$File.Tx[[i]]] %>% GRangesList() %>% unlist() %>% reduce()
})

## find the peaks that are only accessible in a treatment condidtion
treat_only_peaks = mapply(setdiff, treat_peaks, naive_peaks)
names(treat_only_peaks) = df_ifn$Label.Tx
lengths(treat_only_peaks)

## find the cell type peaks constant across treatment and Naive conditions
shared_peaks = mapply(intersect, treat_peaks, naive_peaks)
names(shared_peaks) = df_ifn$Label.Shared
lengths(shared_peaks)

## sanity checks, should be 0's, 0's, lots
mapply(findOverlaps, shared_peaks, treat_only_peaks)
mapply(findOverlaps, naive_peaks, treat_only_peaks)
mapply(findOverlaps, treat_peaks, treat_only_peaks)

###############################################
## 3) export the peaks to make LDSC annotations
treat_only_peaks_fn = here(DATADIR, 'bed', 
                           paste0(names(treat_only_peaks), '.bed.gz'))
mapply(export, treat_only_peaks, treat_only_peaks_fn)

shared_peaks_fn = here(DATADIR, 'bed', paste0(names(shared_peaks), '.bed.gz'))
mapply(export, shared_peaks, shared_peaks_fn)



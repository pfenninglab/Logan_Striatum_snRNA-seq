## packages for data table processing 
library(tidyverse)
library(data.table)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/compare_to_neural_activity_mouse'
dir.create(here(DATADIR, 'rdas'), recursive = T, showWarnings = F)

main_types = c('Neuron', 'Astrocytes', 'Microglia','Oligos',  'Oligos_Pre')
rename = c('Astro' = 'Astrocytes', 'Exc' = 'Neuron', 'Int' = 'Neuron', 
           'Micro' = 'Microglia', 'Olig' = 'Oligos', 'OPC' = 'Oligos_Pre')

mm_to_hg_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv') %>% 
  rename_with(make.names) %>% 
  dplyr::rename('gene' = 'Gene.name','Ensembl' = 'Mouse.gene.stable.ID' )
  

################################################################################
## 1) read in the Blanchard2022_41586_2022_5439_MOESM17_ESM.xlsx APOE-4 study DEGs
hrvatin_table_fn = here(DATADIR, 'tables', 'Hrvatin2017_41593_2017_29_MOESM5_ESM_S3.xlsx')
direction_recode =  setNames(rep(c('up','down'), times = 2), letters[1:4])
time_recode = setNames(rep(c('Early', 'Late'), each = 2), letters[1:4])

## the first sheet has the direction of the foldchange in early or late genes
df_dir = hrvatin_table_fn %>% 
  readxl::read_xlsx(sheet = 'Activity Regulated Genes', skip = 6) %>% 
  dplyr::rename('Mouse.gene.name' = 'Gene Name') %>% 
  pivot_longer(-c(Mouse.gene.name, Ensembl:secreted), names_to = 'celltype', 
               values_drop_na = T, values_to ='direction') %>% 
  dplyr::select(-c(TF, secreted)) %>% 
  mutate(group = time_recode[direction], 
         direction = direction_recode[direction], 
         celltype = celltype %>% ss('_') %>% str_replace('L[1-6]*$', ''),
         celltype = rename[celltype]) %>% distinct()

table(df_dir$celltype)

## the next 2 sheets are the q-values for the DEGs in early or late genes
q_value_sheets = c('Early' = '0h vs 1h q values','Late' = '0h vs 4h q values')
df_qval = lapply(q_value_sheets,  readxl::read_xlsx, skip = 2, path = hrvatin_table_fn) %>% 
  rbindlist(idcol = 'group') %>% dplyr::rename('Mouse.gene.name' = `...1`) %>% 
  pivot_longer(-c(Mouse.gene.name, group), names_to = 'celltype', 
               values_drop_na = T, values_to = 'qval') %>% 
  mutate(celltype = celltype %>% ss('_') %>% str_replace('L[1-6]*$', ''),
         celltype = rename[celltype], 
         qval = as.numeric(qval)) %>% distinct() %>% filter(!is.na(qval))
  
## combine the 2 tables and add the mouse to human genes
df = inner_join(df_dir, df_qval)%>% inner_join(mm_to_hg_genes) %>% 
  filter(!is.na(celltype))
saveRDS(df, here(DATADIR, 'rdas', 'Hrvatin2017_41593_2017_29_MOESM5_ESM_S3.rds'))


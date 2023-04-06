## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(rcartocolor)
library(ggsci)

## stats
library(SeuratDisk)
library(tidymodels)
library(broom)
library(broom.mixed)
library(lme4)
library(lmerTest)
library(swfdr)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

DATADIR= 'data/tidy_data/AUCell_gene_set_activities'
PLOTDIR= 'figures/exploratory/AUCell_gene_set_activities'

###################################
# 0) pre-set plotting aesthetics
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

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

in2mm<-25.4
my_theme = theme_classic(base_size = 6)
alpha = 0.05

#####################################################################
# 1)  load in the filtered Seurat object and grab the cell metadata 
meta_df = readRDS(here('data/tidy_data/AUCell_gene_set_activities', 
                       'rdas/BU_OUD_Striatum_refined_all_SeuratObj_N22.metadata.rds')) %>% 
  dplyr::select(-contains(c('integrated', 'index', 'Index'))) %>% 
  rownames_to_column('Cell')

pd = readRDS('data/tidy_data/tables/LR_RM_OUD_snRNAseq_SampleInfo.rds')
pd2 = pd %>% dplyr::select(Case, ID, Region, Sex, DSM.IV.OUD, PMI, Age) %>% distinct()

############################################################
# 2) read in the import list of inferred regulon target genes
aucell_df = here(DATADIR, 'regulon/OUD_Striatum_refined_N22.aucell.tsv') %>% 
  fread() %>% rename_with(~str_replace(.x, '\\(\\+\\)', '')) %>% 
  pivot_longer(cols = -Cell, names_to = 'TF', values_to = 'AUCell') 

aucell_df = aucell_df %>% inner_join(x= meta_df, y = .) %>% 
  filter(!is.na(AUCell))

table(aucell_df$celltype3, aucell_df$TF)

## compute pseudo-bulk average of AUCell scores by cell type
aucell_pb = aucell_df %>% group_by(ID, Case, Region, TF, celltype3) %>% 
  summarize(AUCell.sd = sd(AUCell, na.rm = T), 
            AUCell = mean(AUCell, na.rm = T)) %>% 
  inner_join(x = pd2, y = .) %>% 
  inner_join(AveCelltypeExpr) %>% 
  filter(!is.na(AUCell.sd)) %>% 
  filter(celltype3 != 'Mural') %>% 
  filter(AUCell.sd >0) %>% 
  group_by(TF, celltype3) %>% 
  ## filter regulons w/o any AUC
  filter(n() >= 10) %>% filter(length(unique(DSM.IV.OUD))>1) %>% 
  filter(length(unique(Sex))>1) %>%ungroup()
  
table(aucell_pb$TF)
table(aucell_pb$celltype3)

## save the pseudo bulk TF regulon output 
save_pb_fn = here(PLOTDIR, 'rdas', 'OUD_Striatum_pyscenicl_TF_modules.pseudoBulkByCelltype.rds')
saveRDS(aucell_pb, save_pb_fn)

##########################################
## 3) mixed effect statistics by cell type
stats_group = aucell_pb %>% 
  nest(data = -c(TF, celltype3)) %>% 
  mutate(
    test = map(data, ~ lm(AUCell ~ DSM.IV.OUD + Age + Sex + PMI + Region, 
                          data = .x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  select(-data, -test) %>% 
  filter(grepl('DSM.IV.OUD', term)) %>% 
  mutate(p.bonferroni = p.adjust(p.value, 'bonferroni'), 
         p.adj = p.adjust(p.value, 'fdr')) %>% 
  arrange(p.value) %>% 
  as.data.frame()

head(stats_group)
stats_group %>% filter(p.adj< alpha)
stats_group %>% filter(p.adj< alpha) %>% count(celltype3)
stats_group %>% filter(p.adj< alpha) %>% count(TF)

stats_group %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_pyscenic_differential_TF_modules.byCelltype.xlsx'))

####################################################################
## 4) compute pseudo-bulk average of AUCell scores by subject
aucell_pb2 = aucell_df %>% group_by(ID, TF) %>% 
  summarize(AUCell.sd = sd(AUCell, na.rm = T), 
            AUCell = mean(AUCell, na.rm = T)) %>% 
  inner_join(x = pd2, y = ., multiple = "all") %>% 
  filter(!is.na(AUCell.sd)) %>% 
  filter(AUCell.sd >0) %>% 
  group_by(TF) %>% 
  ## filter regulons w/o any AUC
  filter(n() >= 10) %>% filter(length(unique(DSM.IV.OUD))>1) %>% 
  filter(length(unique(Sex))>1) %>%ungroup()

table(aucell_pb2$TF)

##  mixed effect statistics by person
stats_group2 = aucell_pb2 %>% 
  nest(data = -c(TF)) %>% 
  mutate(
    test = map(data, ~ lm(AUCell ~ DSM.IV.OUD + Age + Sex + PMI + Region, 
                          data = .x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  select(-data, -test) %>% 
  filter(grepl('DSM.IV.OUD', term)) %>% 
  mutate(p.bonferroni = p.adjust(p.value, 'bonferroni'), 
         p.adj = p.adjust(p.value, 'fdr')) %>% 
  arrange(p.value) %>% 
  as.data.frame()

head(stats_group2)

stats_group2 %>% filter(p.adj< alpha)
stats_group2 %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_pyscenic_differential_TF_modules.bySample.xlsx'))



  
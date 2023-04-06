library(tidyverse)
library(data.table)
library(ggpubr)
library(metafor)
library(tidymeta)
library(broom)
library(qvalue)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/ldsc_interferon_response'
PLOTDIR='figures/exploratory/ldsc_interferon_response'
dir.create(here(DATADIR, 'bed'), showWarnings = F, recursive = T)
dir.create(here(DATADIR, 'rdas'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR, 'plots'), showWarnings = F, recursive = T)
dir.create(here(PLOTDIR, 'tables'), showWarnings = F, recursive = T)

#########################################
## 1) read in the annotation comparisons
# df_ifn = readRDS(here(DATADIR, 'rdas', 'interferon_macrophages_monocyte_annotations.rds')) %>% 
#   dplyr::select(-starts_with('File'), -Group) %>% 
#   pivot_longer(cols = c(Label.Tx, Label.Shared), names_to = 'Treatment', values_to = 'Label') %>% 
#   mutate(Treatment = ifelse(Treatment == 'Label.Tx', paste0(Treatment, '_Unique'), 'Naive_Shared')) %>% 
#   dplyr::rename('Name' = 'Label')


df_ifn = readRDS(here(DATADIR, 'rdas', 'interferon_macrophages_monocyte_annotations.rds')) %>% 
  dplyr::select(-starts_with('File'), -Group) %>% 
  mutate(Label.Tx = str_replace(Label.Tx, '_Unique', '')) %>% 
  dplyr::rename('Name' = 'Label.Tx')


#################################
## 2) read in the list of GWAS
keep_group = c('Degen', 'Other', 'SU', 'SUD', 'Neuro')
pheno = file.path('/projects/pfenninggroup/machineLearningForComputationalBiology', 
                  'gwasEnrichments/gwas_list_sumstats.tsv') %>% 
  read_tsv() %>% filter(group %in% keep_group) %>% filter(!is.na(match)) %>% 
  filter(subgroup != 'Pain') %>% dplyr::select(-file)
table(pheno$group)

################################################
## 3) find the LDSC tau coefficient enrichments
alpha = 0.05
enrich_fn = list.files(here(DATADIR, 'enrichments_naive_bg'),pattern = 'cell_type_results.txt',
                       full.names = T)
names(enrich_fn) = enrich_fn %>% basename() %>% ss('\\.', 2)
enrich_df = lapply(enrich_fn, fread) %>% rbindlist(idcol = 'match') %>% 
  inner_join(x = pheno) %>% inner_join(x= df_ifn) %>% 
  mutate(Coefficient_FDR = qvalue(Coefficient_P_value, fdr.level = alpha)$qvalues) %>% 
  arrange(Coefficient_P_value)

enrich_df %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'gwas_enrichment_interferon_ldsc_tau_coefficient.xlsx'))


enrich_df2 = enrich_df %>% filter(Source != 'Pfenning lab')
enrich_meta = enrich_df2 %>% nest(data = -c(match, Treatment)) %>% 
  mutate(meta = map(data, ~ rma(method = 'EE', data = .x, 
                                yi = Coefficient, sei = Coefficient_std_error)),
         tidied = map(meta, tidy)) %>% 
    unnest(c(tidied)) %>% dplyr::select(-c(data)) %>% 
    mutate(Coefficient_FDR = qvalue(p.value, fdr.level = alpha)$qvalues) %>% 
    dplyr::rename('Coefficient' = 'estimate', 
                  'Coefficient_std_error' = 'std.error') %>% arrange(Coefficient_FDR)

enrich_meta %>% 
  writexl::write_xlsx(here(PLOTDIR, 'tables', 'gwas_meta_interferon_ldsc_tau_coefficient.xlsx'))


################################################
## 4) plot the enrichment
plot_fn = here(PLOTDIR, 'plots', 'interferon_ldsc_tau_SU.pdf')
pdf(plot_fn, height = 6, width = 8)
ggplot(enrich_df %>% filter(group %in% c('SU')), 
       aes(x = Source, y = Coefficient, fill = Treatment)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=Coefficient-Coefficient_std_error, 
                    ymax=Coefficient+Coefficient_std_error), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(subgroup~match %>% str_replace('-', ' ') %>% str_wrap(width = 6)) + coord_flip() + 
  theme_bw(base_size = 6) + 
  theme(legend.position = 'bottom')
dev.off()


plot_fn = here(PLOTDIR, 'plots', 'interferon_ldsc_tau_SUD.pdf')
pdf(plot_fn, height = 6, width = 8)
ggplot(enrich_df %>% filter(group %in% c('SUD')), 
       aes(x = Source, y = Coefficient, fill = Treatment)) +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  geom_errorbar(aes(ymin=Coefficient-Coefficient_std_error, 
                    ymax=Coefficient+Coefficient_std_error), width=.2,
                position=position_dodge(.9)) +
  facet_wrap(subgroup~match %>% str_replace('-', ' ') %>% str_wrap(width = 6)) + coord_flip() + 
  theme_bw(base_size = 6) + 
  theme(legend.position = 'bottom')
dev.off()




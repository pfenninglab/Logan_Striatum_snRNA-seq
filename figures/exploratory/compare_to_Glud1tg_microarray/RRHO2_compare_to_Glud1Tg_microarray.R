ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(data.table)
library(tidyverse)
library(RRHO2)
library(here)
library(broom)

DATADIR='data/tidy_data/compare_to_Glud1tg_microarray'

## make for this subdirs
PLOTDIR='figures/exploratory/compare_to_Glud1tg_microarray'
here(PLOTDIR, c('plots/RRHO2', 'tables', 'rdas')) %>% 
  sapply(dir.create, showWarnings = F)

##################################################################
# 1) grab the OUD DEGs split by the DEGs to plot in neurons or glia
res = here('data/tidy_data/differential_expression_analysis', 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% 
  readRDS()
oud_deg_list = res %>% 
  lapply(function(x){
    x %>% mutate(value = sign(logFC) * -log10(P.Value)) %>% 
      dplyr::select(gene, value) %>% filter(!is.na(value) & !duplicated(gene))
  }) %>% lapply(as.data.frame)

######################################################
# 1) load in the DEGs from Glud1Tg microarray analyses 
deg_df = readRDS(here(DATADIR, 'rdas', 'Glud1Tg_microarray_reanalysis_GSE48911_GSE11419.rds')) %>% mutate(value = sign(logFC) * -log10(P.Value)) %>% dplyr::rename('gene' = 'Gene.name') %>% 
  filter(!duplicated(gene))
deg_list = rep(list(deg_df), length(oud_deg_list))
names(deg_list) =names(oud_deg_list)

############################
# 3) run the RRHO2 objects
ind_list = set_names(names(deg_list))
deg_list = lapply(ind_list, function(ind){
  shared = intersect(deg_list[[ind]]$gene, oud_deg_list[[ind]]$gene)
  deg_list[[ind]] %>% filter(gene %in% shared) %>% 
  dplyr::select(gene, value)
})
oud_deg_list = lapply(ind_list, function(ind){
  shared = intersect(deg_list[[ind]]$gene, oud_deg_list[[ind]]$gene)
  oud_deg_list[[ind]] %>% filter(gene %in% shared)
})

## create the Rank-rank hypermetric overlap objects
RRHO_list = parallel::mclapply(ind_list, function(ind){
  RRHO2_initialize( list1 = oud_deg_list[[ind]], list2 = deg_list[[ind]],
  labels = c("OUD DEG", ind), boundary = 0.05, log10.ind=TRUE)
}, mc.cores = 12)

RRHO_list %>% saveRDS(here(PLOTDIR,'rdas', paste0('compare_to_Glud1tg_microarray_rrho.rds')))

out = lapply(ind_list, function(ind){
  df = inner_join(x = oud_deg_list[[ind]] %>% rename('oud' = 'value'), 
                  y =  deg_list[[ind]] %>% rename('glud1' = 'value'))
  cor.test( df$oud , df$glud1, method = 'spearman') %>% tidy()
}) %>% rbindlist(idcol = 'celltype') %>% arrange(p.value) %>% 
  mutate(fdr = p.adjust(p.value))

############################################################
# 4) plot the RRHO2 grouped all together with one scale
in2mm<-25.4
my_theme = theme_classic(base_size = 6)
RRHO_list = RRHO_list[!is.na(RRHO_list %>% lapply('[', 'hypermat'))]
maximum <- RRHO_list %>% lapply('[', 'hypermat') %>% unlist() %>% max(na.rm=TRUE)
minimum <- RRHO_list %>% lapply('[', 'hypermat') %>% unlist() %>% min(na.rm=TRUE)

plot_fn = here(PLOTDIR,'plots', paste0('compare_to_Glud1tg_microarray_rrho.pdf'))
pdf(plot_fn, height = 100/in2mm, width = 100/in2mm)
for(ind in names(RRHO_list)){
  par(mar = rep(1.65, 4))
  RRHO2_heatmap(RRHO_list[[ind]], maximum=maximum, minimum=minimum,
                labels = c( paste(ss(ind, '#'), 'OUD DEG'), 'Glud1Tg Mouse DEG'))
}
dev.off()

############################################################
# 5) do the gene set overlap w/ a simple rank correlation
out %>% dplyr::select(-c(statistic, method, alternative)) %>% 
  dplyr:: rename('Spearmans Rho'= 'estimate')

out %>% saveRDS(here(PLOTDIR,'rdas', paste0('compare_to_Glud1tg_microarray_spearman.rds')))
out %>% writexl::write_xlsx(here(PLOTDIR,'tables', paste0('compare_to_Glud1tg_microarray_spearman.xlsx')))



ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(data.table)
library(tidyverse)
library(RRHO2)
library(here)

DATADIR='data/tidy_data/compare_to_Alzheimers_studies'

## make for this subdirs
PLOTDIR='figures/explanatory/figureNN_compare_alzheimers'
here(PLOTDIR, c('plots/RRHO2', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

################################################################
# 1) load in the DEGs from AD studies grouped by neuron and glia

## Blanchard et al. DEGs from snRNA or sorted/cultured Oligos
main_types = c('Neuron', 'Astrocytes', 'Microglia','Oligos',  'Oligos_Pre')
ad_deg_list1 = readRDS(here(DATADIR, 'rdas', 'Blanchard2022_APOE_allele_celltype_DEG_list.rds')) %>% 
  lapply(function(x){ x %>% mutate(value = sign(logFC) * -log10(P.Value))})

ad_deg_list2 = readRDS(here(DATADIR, 'rdas', 'Blanchard2022_APOE_PMB_Oligo_DEG_list.rds')) %>% 
  lapply(function(x){ x %>% mutate(value = sign(logFC) * -log10(pval)) })

## Mathys et al. DEGs
ad_deg_list3 = readRDS(here(DATADIR, 'rdas', 'Mathys2018_41586_2019_1195_MOESM4_ESM.rds')) %>% 
  lapply(function(x){
    x %>% mutate(value = sign(as.numeric(IndModel.FC)) * -log10(as.numeric(IndModel.adj.pvals))) })

## combine the lists together
ad_deg_list = c(ad_deg_list1, ad_deg_list2, ad_deg_list3) %>% 
  lapply(as.data.frame) %>% lapply(function(x){
    x  %>% dplyr::select(gene, value) %>% filter(!is.na(value) & !duplicated(gene))
  })

ad_deg_list = ad_deg_list[sort(names(ad_deg_list))]
ad_deg_list %>% sapply(function(x) x %>% filter(gene =='MTRNR2L8')) %>% t()
ad_deg_list %>% sapply(function(x) x %>% filter(gene =='APOE')) %>% t()

##################################################################
# 2) grab the OUD DEGs split by the DEGs to plot in neurons or glia
res = here('data/tidy_data/differential_expression_analysis', 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% 
  readRDS()
oud_deg_list = res[main_types]%>% 
  lapply(function(x){
    x %>% mutate(value = sign(logFC) * -log10(P.Value)) %>% 
      dplyr::select(gene, value) %>% filter(!is.na(value) & !duplicated(gene))
  }) %>% lapply(as.data.frame)

indOUD2AD = match(names(ad_deg_list) %>% ss('#'), names(oud_deg_list))
oud_deg_list = oud_deg_list[indOUD2AD]
names(oud_deg_list) = names(ad_deg_list)

############################
# 3) run the RRHO2 objects
ind_list = set_names(names(ad_deg_list))
ad_deg_list = lapply(ind_list, function(ind){
  shared = intersect(ad_deg_list[[ind]]$gene, oud_deg_list[[ind]]$gene)
  ad_deg_list[[ind]] %>% filter(gene %in% shared)
})
oud_deg_list = lapply(ind_list, function(ind){
  shared = intersect(ad_deg_list[[ind]]$gene, oud_deg_list[[ind]]$gene)
  oud_deg_list[[ind]] %>% filter(gene %in% shared)
})

## create the Rank-rank hypermetric overlap objects
RRHO_list = parallel::mclapply(ind_list, function(ind){
  RRHO2_initialize( list1 = oud_deg_list[[ind]], list2 = ad_deg_list[[ind]],
  labels = c("OUD DEG", ind), boundary = 0.05, log10.ind=TRUE)
}, mc.cores = 12)

RRHO_list %>% saveRDS(here(PLOTDIR,'rdas', paste0('compare_alzheimers_rrho.rds')))

############################################################
# 4) plot the RRHO2 grouped all together with one scale
in2mm<-25.4
my_theme = theme_classic(base_size = 6)
RRHO_list = RRHO_list[!is.na(RRHO_list %>% lapply('[', 'hypermat'))]
maximum <- RRHO_list %>% lapply('[', 'hypermat') %>% unlist() %>% max(na.rm=TRUE)
minimum <- RRHO_list %>% lapply('[', 'hypermat') %>% unlist() %>% min(na.rm=TRUE)

plot_fn = here(PLOTDIR,'plots', paste0('figureNN_compare_alzheimers_rrho.pdf'))
pdf(plot_fn, height = 100/in2mm, width = 100/in2mm)
for(ind in names(RRHO_list)){
  par(mar = rep(1.65, 4))
  RRHO2_heatmap(RRHO_list[[ind]], maximum=maximum, minimum=minimum,
                labels = c( paste(ss(ind, '#'), 'OUD DEG'), gsub('#|_', ' ', ind)))
}
dev.off()


############################################################
# 5) plot the RRHO2 grouped all together with separate scales
ind_list2 = split(names(RRHO_list), grepl('APOE', names(RRHO_list)))
names(ind_list2) = c('Mathys2018', 'Blanchard2022')
for(lab in names(ind_list2)){
  ## find the range per study
  maximum <- RRHO_list[ind_list2[[lab]]] %>% 
    lapply('[', 'hypermat') %>% unlist() %>% max(na.rm=TRUE)
  minimum <-  RRHO_list[ind_list2[[lab]]] %>% 
    lapply('[', 'hypermat') %>% unlist() %>% min(na.rm=TRUE)
  
  ## plot the data split by AD study
  plot_fn = here(PLOTDIR,'plots', paste0('figureNN_compare_alzheimers_rrho.',lab,'.pdf'))
  pdf(plot_fn, height = 100/in2mm, width = 100/in2mm)
  for(ind in ind_list2[[lab]]){
    plot_fn = here(PLOTDIR,'plots', 'RRHO2', paste0('figureNN_compare_alzheimers_rrho.',ind,'.pdf'))
    pdf(plot_fn, height = 100/in2mm, width = 100/in2mm)
    par(mar = rep(1.65, 4))
    RRHO2_heatmap(RRHO_list[[ind]], maximum=maximum, minimum=minimum,
                  labels = c(paste(ss(ind, '#'), 'OUD DEG'), gsub('#|_', ' ', ind)))
    dev.off()
    par(mar = rep(1.65, 4))
    RRHO2_heatmap(RRHO_list[[ind]], maximum=maximum, minimum=minimum,
                    labels = c(paste(ss(ind, '#'), 'OUD DEG'), gsub('#|_', ' ', ind)))
  }
  dev.off()
}
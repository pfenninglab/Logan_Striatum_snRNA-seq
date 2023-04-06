ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(data.table)
library(tidyverse)
library(RRHO2)
library(here)

DATADIR='/projects/pfenninggroup/singleCell/McLean_chronic_opioid_monkey_snRNA-seq/data/tidy_data/differential_expression'

## make for this subdirs
PLOTDIR='figures/exploratory/compare_McLean_monkey_snRNA-seq'
here(PLOTDIR, c('plots/RRHO2', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F, recursive = T)

###################################################################
# 1) load in the DEGs from monkey snRNA-seq NAc exposed to opioids
main_types = c('All','Neuron', 'Glia', 'Interneuron','Astrocytes', 'Microglia','Oligos',  'Oligos_Pre')

celltype_mappings = 
  monkey_deg_list = here(DATADIR, 'rdas', 'Monkey_voom_limma_bigModelSVA_N16.celltype.rds') %>% 
  readRDS() %>% 
  lapply(function(x){
    x %>% mutate(value = sign(logFC) * -log10(P.Value)) %>% 
      dplyr::select(gene, value) %>% filter(!is.na(value) & !duplicated(gene))
  })

monkey_deg_list %>% sapply(function(x) x %>% filter(gene =='APOE')) %>% t()

##################################################################
# 2) grab the OUD DEGs split by the DEGs to plot in neurons or glia
res = here('data/tidy_data/differential_expression_analysis', 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% 
  readRDS()
oud_deg_list = res %>% 
  lapply(function(x){
    x %>% mutate(value = sign(logFC) * -log10(P.Value)) %>% 
      dplyr::select(gene, value) %>% filter(!is.na(value) & !duplicated(gene))
  }) %>% lapply(as.data.frame)


####################################################################################
# 3) align the cell type DEGs from monkey with cell type DEGs from human w/ first 3 characters
ind_list = paste(rep(names(monkey_deg_list), each = length(oud_deg_list)), 
                rep(names(oud_deg_list), times = length(monkey_deg_list)), sep = '#')
names(ind_list) = ind_list

monkey_deg_list2 = lapply(ind_list, function(ind){
  i1 = ss(ind, '#', 1); i2 = ss(ind, '#', 2)
  shared = intersect(monkey_deg_list[[i1]]$gene, oud_deg_list[[i2]]$gene)
  monkey_deg_list[[i1]] %>% filter(gene %in% shared) %>% 
    filter(!duplicated(gene), !is.na(gene))
})
oud_deg_list2 = lapply(ind_list, function(ind){
  i1 = ss(ind, '#', 1); i2 = ss(ind, '#', 2)
  shared = intersect(monkey_deg_list[[i1]]$gene, oud_deg_list[[i2]]$gene)
  oud_deg_list[[i2]] %>% filter(gene %in% shared) %>% 
    filter(!duplicated(gene), !is.na(gene))
})

## create the Rank-rank hypermetric overlap objects
RRHO_list = parallel::mclapply(ind_list, function(ind){
  RRHO2_initialize( list1 = oud_deg_list2[[ind]], list2 = monkey_deg_list2[[ind]],
  labels = c("OUD DEG", ind), boundary = 0.05, log10.ind=TRUE)
}, mc.cores = 12)

RRHO_list %>% saveRDS(here(PLOTDIR,'rdas', paste0('compare_McLean_monkey_snRNA-seq_rrho.rds')))

############################################################
# 4) plot the RRHO2 grouped all together with one scale
in2mm<-25.4
my_theme = theme_classic(base_size = 6)

## pare down RRHO's to the only matching cell type comparisons
ind_list2 = ind_list[ss(ind_list, '#', 1) %>% str_sub(1,3) ==
                       ss(ind_list, '#', 2) %>% str_sub(1,3)]

# drop these cell type combinations, faulty matches
ind_list2 = ind_list2[!ind_list2 %in% c("D1.ICj#D1.D2.Hybrid", "D1.ICj#D1.Matrix",
                                        "D1.ICj#D1.Striosome", "D1.NUDAP#D1.Matrix",
                                        "D1.NUDAP#D1.Striosome", "D1.Shell.OT#D1.D2.Hybrid",
                                        "Oligos_Pre#Oligos", "Oligos#Oligos_Pre")]

RRHO_list = RRHO_list[!is.na(RRHO_list %>% lapply('[', 'hypermat'))]
maximum <- RRHO_list[ind_list2] %>% lapply('[', 'hypermat') %>% unlist() %>% max(na.rm=TRUE)
minimum <- RRHO_list[ind_list2] %>% lapply('[', 'hypermat') %>% unlist() %>% min(na.rm=TRUE)


plot_fn = here(PLOTDIR,'plots', paste0('figureNN_compare_McLean_monkey_snRNA-seq_rrho.pdf'))
pdf(plot_fn, height = 100/in2mm, width = 100/in2mm)
for(ind in ind_list2){
  par(mar = rep(1.65, 4))
  i1 = ss(ind, '#', 1); i2 = ss(ind, '#', 2)
  RRHO2_heatmap(RRHO_list[[ind]], maximum=maximum, minimum=minimum,
                labels = c( paste(i2, 'hsOUD DEG'), paste(i1, 'rmOUD DEG')))
}
dev.off()



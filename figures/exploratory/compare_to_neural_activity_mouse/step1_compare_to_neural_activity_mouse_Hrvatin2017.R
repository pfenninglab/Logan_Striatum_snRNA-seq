## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

library(data.table)
library(fgsea)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


## make for this subdirs
DATADIR='data/tidy_data/compare_to_neural_activity_mouse'
PLOTDIR='figures/exploratory/compare_to_neural_activity_mouse'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

main_types = c('Neuron', 'Astrocytes', 'Microglia','Oligos',  'Oligos_Pre')
rename = c('Astro' = 'Astrocytes', 'Exc' = 'Neuron', 'Int' = 'Neuron', 
           'Micro' = 'Microglia', 'Olig' = 'Oligos', 'OPC' = 'Oligos_Pre')

##############################
# 1) read in the DEG lists
res.celltype = here('data/tidy_data/differential_expression_analysis/rdas', 
                    'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% readRDS()
names(res.celltype) = paste0(names(res.celltype), '#All')

deg_rank_list = lapply(c(res.celltype), function(deg){
  deg %>% mutate(tmp = -log10(P.Value) * sign(logFC)) %>% 
    dplyr::select(gene, tmp) %>% arrange(tmp) %>% deframe()
})

## read in the gene sets from the Hrvatin, Hochbaum 2017 et al. paper
hrvatin_list = readRDS(here(DATADIR, 'rdas', 'Hrvatin2017_41593_2017_29_MOESM5_ESM_S3.rds')) %>% 
  mutate(split = paste(celltype, group, direction, sep = '#')) %>% 
  split(x = .$gene, f = .$split) %>% lapply(unique)
lengths(hrvatin_list)
names(hrvatin_list)[lengths(hrvatin_list)>=5]

## conduct the GSEA analyses
gsea_list = lapply(deg_rank_list, fgsea, pathways = hrvatin_list,
                   minSize=5, ## minimum gene set size
                   maxSize=400) ## maximum gene set size

gsea_df = gsea_list %>% rbindlist(idcol = 'group') %>% arrange(pval) %>% 
  mutate(
    celltype = group %>% ss('#', 1), 
    celltype2 = pathway %>% ss('#', 1),
    direction = pathway %>% ss('#', 2),
    OUD.v.CTL.in = group %>% ss('#', 2), 
    filt = case_when(grepl("Neu|D1|D2|Int", celltype) & celltype2== 'Neuron' ~T,
                     celltype == celltype2 ~ T, T ~ F),
    leadingEdge = map_chr(leadingEdge, paste, collapse = ',')) %>% 
  filter(filt) %>% dplyr::select(-filt) %>% 
  mutate(padj = p.adjust(pval, 'fdr')) %>% 
  dplyr::select(-group) %>% relocate(celltype, OUD.v.CTL.in, .before= everything())

out_fn = here(PLOTDIR, 'tables', 'GSEA_enrichment_hrvatin_et_al_mouse_v1c_light_exposed.xlsx')
gsea_df %>% writexl::write_xlsx(out_fn)


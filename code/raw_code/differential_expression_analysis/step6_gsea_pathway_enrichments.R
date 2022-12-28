## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

library(data.table)
library(fgsea)
library(swfdr)
library(msigdbr)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)


## make for this subdirs
PLOTDIR='figures/exploratory/Piechota_et_al_human_striatum_on_drugs'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)


############################################################
# 1) read in the DEG lists per comparison of OUD vs. Control
res.celltype = here('data/tidy_data/differential_expression_analysis/rdas', 
                    'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% readRDS()
names(res.celltype) = paste0(names(res.celltype), '#All')

res.OUDwinSex = here('data/tidy_data/differential_expression_analysis/rdas', 
                     'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds') %>% 
  readRDS() %>% rbindlist() %>% dplyr::select(-contains('dir')) %>% 
  pivot_longer(cols = -c(celltype, gene), names_sep = '_Sex', names_to = c('name', 'Sex')) %>% 
  pivot_wider(id_cols =c(celltype, gene, Sex)) %>% 
  mutate(group = paste0(celltype, '#Sex', Sex)) %>% 
  split(f = .$group)

res.OUDwinRegion = here('data/tidy_data/differential_expression_analysis/rdas', 
                        'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinRegion.rds') %>% 
  readRDS() %>% rbindlist() %>% dplyr::select(-contains('dir')) %>% 
  pivot_longer(cols = -c(celltype, gene), names_sep = '_Region', names_to = c('name', 'Region')) %>% 
  pivot_wider(id_cols =c(celltype, gene, Region)) %>% 
  mutate(group = paste0(celltype, '#Region', Region)) %>% 
  split(f = .$group)

deg_rank_list1 = lapply(c(res.celltype, res.OUDwinSex, res.OUDwinRegion), function(deg){
  deg %>% mutate(tmp = -log10(P.Value) * sign(logFC)) %>% 
    dplyr::select(gene, tmp) %>% arrange(tmp) %>% deframe()
})

## read in the DEG columns again using the interaction columns
res.OUDwinSex = here('data/tidy_data/differential_expression_analysis/rdas', 
                     'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds') %>% 
  readRDS() %>% rbindlist() %>%
  mutate(group = paste0(celltype, '#SexInteraction')) %>% 
  split(f = .$group)

res.OUDwinRegion = here('data/tidy_data/differential_expression_analysis/rdas', 
                        'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinRegion.rds') %>% 
  readRDS() %>% rbindlist() %>%
  mutate(group = paste0(celltype, '#RegionInteraction')) %>% 
  split(f = .$group)

deg_rank_list2 = lapply(c(res.OUDwinSex, res.OUDwinRegion), function(deg){
  deg %>% dplyr::select(gene, dir_difference) %>% arrange(dir_difference) %>% deframe()
})

## combine the 2 ranked DEG by cell type comparison lists
deg_rank_list = c(deg_rank_list1, deg_rank_list2)



#############################################
## 2) get gene ontologies, use human genes
## grab the H, Hallmark set of gene pathways
## grab the C2, which is the curated canonical pathway sets
## grab the C5, which is the Gene Ontology sets
## https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
pathways_df =  bind_rows(msigdbr("human", category="H"), 
                         msigdbr("human", category="C2"), 
                         msigdbr("human", category="C5"))

## get the SynGO gene ontologies, use human genes
syngo_df = readRDS(here('data/tidy_data/SynGO_bulk_download_release_20210225',
                           'rdas', 'syngo_annotations.rds'))
pathways_df = rbindlist(list(pathways_df, syngo_df), fill = T) 

## reshape/label for pathway naming purposes
pathways <-pathways_df %>% 
  mutate(
    gs_subcat = ifelse(is.na(gs_subcat) | gs_subcat == '', gs_cat, gs_subcat),
    gs_name = paste(gs_subcat, gs_name, sep ='#')) %>% 
  split(x = .$gene_symbol, f = .$gs_name)


## exclude the really really big gene sets
lengths(pathways) %>% summary()
pathways = pathways[lengths(pathways)<500]
length(pathways) # 21308

table(pathways_df$gs_cat)

pathways_df2 = pathways_df %>% 
  dplyr::select(gs_subcat, gs_name, gs_description) %>% distinct() %>% 
  dplyr::rename('pathway_group' = 'gs_subcat', 'pathway' = 'gs_name',
                'description' = 'gs_description')


## conduct the GSEA analyses
gsea_list = lapply(deg_rank_list, fgsea, pathways = pathways,
                   minSize=15, ## minimum gene set size
                   maxSize=400) ## maximum gene set size

alpha = 0.05
gsea_df = gsea_list %>% rbindlist(idcol = 'group') %>% 
  arrange(pval) %>% filter(!is.na(pval)) %>% 
  mutate(
    MSigDb_Group = ss(pathway, '#', 1), 
    celltype = group %>% ss('#', 1), 
    OUD.v.CTL.in = group %>% ss('#', 2),
    pathway = ss(pathway, '#', 2),
    padj = lm_qvalue(pval, X=size)$q, 
    celltype = group %>% ss('#', 1), 
    leadingEdge = map_chr(leadingEdge, paste, collapse = ',')) %>% 
  inner_join(pathways_df2) %>% filter(padj < alpha) %>% dplyr::select(-group) %>% 
  relocate(OUD.v.CTL.in, celltype, MSigDb_Group, description, .before= everything()) %>% 
  split(f = .$OUD.v.CTL.in)

sapply(gsea_df, nrow)

out_fn = here('data/tidy_data/differential_expression_analysis', 
              'tables/OUD_Striatum_GSEA_enrichment_msigdb_H_C2_C5_SynGO.xlsx')
gsea_df %>% writexl::write_xlsx(out_fn)
gsea_df %>%saveRDS(here('data/tidy_data/differential_expression_analysis', 
     'rdas/OUD_Striatum_GSEA_enrichment_msigdb_H_C2_C5_SynGO.rds'))

## take a look at the enrichments
gsea_df2 = gsea_df %>% rbindlist()
table(gsea_df2$OUD.v.CTL.in, gsea_df2$MSigDb_Group)

## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(data.table)

## stats
library(fgsea)
library(swfdr)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

DATADIR= 'data/tidy_data/AUCell_gene_set_activities'
PLOTDIR= 'figures/exploratory/AUCell_gene_set_activities'

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

############################################################
# 2) read in the import list of inferred regulon target genes
regulon_df = here(DATADIR, 'regulon/OUD_Striatum_refined_N22.regulons.tsv') %>% 
  fread(skip = 1) %>% filter(V1 != 'TF') %>% 
  mutate(TargetGenes = map(TargetGenes, ~ .x %>% str_replace_all("\\[|\\]|'|\\(",'') %>%
                             str_split('\\), ') %>% unlist())) %>% unnest(TargetGenes) %>% 
  mutate(importance = ss(TargetGenes, ', ', 2) %>% str_replace('\\)|\\]', '') %>% as.numeric(), 
         TargetGenes = ss(TargetGenes, ', ', 1) %>% str_replace('\\[', ''))

regulon_df = regulon_df %>% dplyr::select(-c(V2:Context, RankAtMax)) %>% 
  distinct()

regulon_df %>% writexl::write_xlsx(here(PLOTDIR, 'tables', 'OUD_Striatum_pyscenic_deteced_TF_modules.xlsx'))

pathways = split(regulon_df$TargetGenes, regulon_df$V1)
summary(lengths(pathways)) # regulons have 3 to 2500 genes

## conduct the GSEA analyses on regulons discovered by pyscenic
gsea_list = parallel::mclapply(deg_rank_list, fgsea, pathways = pathways,
                   minSize=5, ## minimum gene set size
                   maxSize=2500, mc.cores = 24) ## maximum gene set size, 


alpha = 0.05
gsea_df = gsea_list %>% rbindlist(idcol = 'group') %>% 
  arrange(pval) %>% filter(!is.na(pval)) %>% 
  mutate(
    celltype = group %>% ss('#', 1), 
    OUD.v.CTL.in = group %>% ss('#', 2),
    padj = lm_qvalue(pval, X=size)$q, 
    p.bonferroni = p.adjust(pval, 'bonferroni'), 
    celltype = group %>% ss('#', 1), 
    leadingEdge = map_chr(leadingEdge, paste, collapse = ',')) %>% 
  dplyr::select(-group) %>% 
  relocate(OUD.v.CTL.in, celltype, pathway, p.bonferroni, padj,  .before= everything()) %>% 
  split(f = .$OUD.v.CTL.in)

## save the enrichment w/ all the pathways, significant and otherwise
gsea_df %>%saveRDS(here(PLOTDIR, 'rdas/OUD_Striatum_GSEA_enrichment_pyscenic_regulons.rds'))

## filter out just the significant pathways
gsea_df = gsea_df %>% lapply(filter, padj < alpha)
sapply(gsea_df, nrow)

out_fn = here(PLOTDIR, 'tables/OUD_Striatum_GSEA_enrichment_pyscenic_regulons.xlsx')
gsea_df %>% writexl::write_xlsx(out_fn)
gsea_df %>%saveRDS(here(PLOTDIR,  'rdas/OUD_Striatum_GSEA_enrichment_pyscenic_regulons.rds'))





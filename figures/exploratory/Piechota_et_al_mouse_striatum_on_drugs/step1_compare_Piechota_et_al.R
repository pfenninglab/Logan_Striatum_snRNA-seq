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
PLOTDIR='figures/exploratory/Piechota_et_al_mouse_striatum_on_drugs'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

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

dx_col = setNames(pal_npg(palette = c("nrc"), alpha = 1)(2), c('CTL', 'OUD'))

in2mm<-25.4
my_theme = theme_classic(base_size = 6)


##############################
# 1) read in the DEG lists
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

deg_rank_list = lapply(c(res.celltype, res.OUDwinSex, res.OUDwinRegion), function(deg){
  deg %>% mutate(tmp = -log10(P.Value) * sign(logFC)) %>% 
    dplyr::select(gene, tmp) %>% arrange(tmp) %>% deframe()
})

## read in the 1-1 gene orthologs
hg38_to_mm10_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv', show_col_types = FALSE)%>% 
  rename_with(make.names) %>% dplyr::select(c(Mouse.gene.name, Gene.name))


## read in the gene sets from the Piechota et al. paper
piechota_list = here('data/tidy_data/Piechota_et_al_mouse_striatum_on_drugs/tables',
               "13059_2010_2341_MOESM2_ESM.XLS") %>%
  readxl::read_xls('patterns_gene_expression_values') %>% 
  rename_with(make.names) %>% dplyr::rename('Mouse.gene.name' = 'Gene.Symbol') %>% 
  dplyr::select(Mouse.gene.name:B3_extended) %>% 
  pivot_longer(cols = A:B3_extended,names_to = 'pattern', values_to = 'value' ) %>% 
  filter(!is.na(value) & grepl('ext', pattern)) %>% 
  mutate(pattern = ss(pattern, '_') %>% paste0('Piechota2010.',.)) %>% 
  inner_join(hg38_to_mm10_genes) %>% split(x = .$Gene.name, f = .$pattern)
    

## conduct the GSEA analyses
gsea_list = lapply(deg_rank_list, fgsea, pathways = piechota_list,
                   minSize=15, ## minimum gene set size
                   maxSize=400) ## maximum gene set size

gsea_df = gsea_list %>% rbindlist(idcol = 'group') %>% arrange(pval) %>% 
  mutate(padj = p.adjust(pval, 'fdr'), 
         celltype = group %>% ss('#', 1), 
         OUD.v.CTL.in = group %>% ss('#', 2), 
         leadingEdge = map_chr(leadingEdge, paste, collapse = ',')) %>% 
  dplyr::select(-group) %>% relocate(celltype, OUD.v.CTL.in, .before= everything())

out_fn = here(PLOTDIR, 'tables', 'GSEA_enrichment_Piechota_et_al_mouse_striatum_on_drugs.xlsx')
gsea_df %>% writexl::write_xlsx(out_fn)


## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RRHO2)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'
PLOTDIR='figures/exploratory/differential_expression_analysis/plots'
dir.create(PLOTDIR, showWarnings = F)

##############################################
# 1) load in the Seney et al. OUD DEG tables
tb_fn = list.files('data/tidy_data/Seney2021_OUDvsCTL_bulk_DLPFC_NAc_tables/tables', 
                   pattern = '.xlsx', full.names = T)
names(tb_fn) = basename(tb_fn) %>% ss('\\.', 3)
degList = tb_fn %>% lapply(readxl::read_xlsx) %>% 
  lapply(function(x) {
    x %>% mutate(DDE = -log10(corrected.p) * sign(coefficient)) %>% 
      dplyr::select(gene_symbol, DDE) %>% as.data.frame() %>%
      arrange(gene_symbol)
  })


########################################
# 3) load in the edgeR perCelltype DEG
save_res_fn = here(DATADIR,'rdas', 'BU_Run1_Striatum_edgeRQLFDetRateRes_byCelltype2_N4.rds')
res = readRDS(save_res_fn) 

# filter FDR < 0.05
alpha = 0.05; LFC_threshold = log2(1.2) # 20% difference
res_fil <- lapply(res, function(u) u %>% 
                    dplyr::filter(p_adj < alpha, abs(logFC) > LFC_threshold) %>% 
                    dplyr::arrange(p_adj))

## get number & percent of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
celltype_degList = lapply(res, function(x) {
    x %>% dplyr::rename('gene_symbol' = 'gene') %>%
    mutate(DDE = -log10(p_adj) * sign(logFC)) %>% 
    dplyr::select(gene_symbol, DDE) %>% as.data.frame() %>%
    arrange(gene_symbol)
  })

################################################################################
# 4) plot the rank-rank hypergeometric overlap of the OUD v. CTL DEG profiles
plot_fn = here(PLOTDIR, 'BU_Run1_Striatum_edgeRQLFDetRateRes_N4.RRHO2wSeney2021.pdf')
pdf(plot_fn, height = 6, width = 6)
for(tissue in names(degList)){
for(cell in names(celltype_degList)){
  l1 = degList[[tissue]]; l2 = celltype_degList[[cell]]
  l1 = l1[l1$gene_symbol %in% l2$gene_symbol, ]
  l2 = l2[l2$gene_symbol %in% l1$gene_symbol, ]
  print(paste0('performing RRHO2 for ', tissue, ' and ', cell, '.'))
  RRHO_obj <- RRHO2_initialize(l2, l1, log10.ind=TRUE,
                               labels = c(cell, paste('Seney et al., 2021', tissue)))
  RRHO2_heatmap(RRHO_obj)
}}
dev.off()







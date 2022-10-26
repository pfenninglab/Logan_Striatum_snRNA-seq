## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main differential gene expression package
library(SingleCellExperiment)
library(DelayedArray)
library(HDF5Array)
library(limma)
library(edgeR)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'
h5Dir =here(DATADIR, 'HDF5Array')
dir.create(h5Dir, showWarnings = F)

PLOTDIR='figures/exploratory/differential_expression_analysis/plots'
dir.create(PLOTDIR, showWarnings = F)

##################################################
# 1) load in cell type labels for label transfer
sce = loadHDF5SummarizedExperiment(h5Dir, prefix="OUD_Striatum_refined_all_SeuratObj_N22.h5Seuratassays")

## sore cell type IDs (kids) and sample IDs (sids)
nk <- length(kids <- rlang::set_names(levels(factor(sce$celltype3))))
ns <- length(sids <- rlang::set_names(levels(factor(sce$orig.ident))))

#############################################################################
# 3) use edgeRQLFDetRate to detect per-cell type differential state analyses 
## for ea. cell type, run edgeRQLFDetRate w/ default parameters
## https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
save_res_fn = here(DATADIR,'rdas', 'OUD_Striatum_edgeRQLFDetRateRes_bycelltype3_N22.rds')
res = readRDS(save_res_fn)

# filter FDR < 0.05, |logFC| > 1 & sort by FDR
alpha = 0.05; LFC_threshold = log2(1.2) # 20% difference
res_fil <- lapply(res, function(u) u %>% 
                    dplyr::filter(p_adj < alpha, abs(logFC) > LFC_threshold) %>% 
                    dplyr::arrange(p_adj))

## get number & percent of DE genes per cluster
n_de <- vapply(res_fil, nrow, numeric(1))
cbind(n_de, p_gs = n_de / nrow(sce) * 100)

## get the top 5 genes of each cell type
top_gs <- lapply(res_fil, function(u) u$gene[seq_len(5)])
top_gs = top_gs[sapply(top_gs, function(x) !all(is.na(x)))]
cs_by_k <- split(colnames(sce), sce$celltype3)

## make scatter plot per subject per cell type per DE gene per cell
plot_fn = here(PLOTDIR, 'OUD_Striatum_edgeRQLFDetRateRes_N22.bycelltype3.topDiffGene.pdf')
pdf(plot_fn, height = 4, width = 7.3)
# split cells by cluster
lapply(names(top_gs), function(k) {
  gs <- top_gs[[k]]  # get top gene-hits for cluster k
  gs = gs[!is.na(gs)]
  cs <- cs_by_k[[k]] # subset cells assigned to cluster k
  scater::plotExpression(sce[, cs], features = gs, x = "orig.ident", colour_by = "DSM.IV.OUD", ncol = 5) +
    guides(fill = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    theme_classic() + ggtitle(paste('Celltype:', k)) + 
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1))
})
dev.off()



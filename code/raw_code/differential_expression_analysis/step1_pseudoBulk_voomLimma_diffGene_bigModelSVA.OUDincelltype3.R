## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(rlang)
library(writexl)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## main differential gene expression package
library(SingleCellExperiment)
library(DelayedArray)
library(HDF5Array)
library(Matrix.utils)
library(limma)
library(edgeR)
library(sva)
library(swfdr)

## regress out the surrogate variables
library(jaffelab)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'

#######################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 28 cores
plan("multicore", workers = 28)
options(future.globals.maxSize = 80 * 1024^3)


##################################################
# 1) create or load pseudobulk sce object
save_pseudobulk =here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds')
if(!file.exists(save_pseudobulk)){
  ## load the single cell counts
  h5Dir =here(DATADIR, 'HDF5Array'); dir.create(h5Dir, showWarnings = F)
  sce = loadHDF5SummarizedExperiment(h5Dir, prefix="BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat")
  
  ## merge interneurons again
  sce$celltype3 = ifelse(grepl('Int', sce$celltype3), 'Interneuron',sce$celltype3)
  table(sce$celltype3)
  
  ## aggregate by cluster-sample to create pseudo-bulk count matrix
  colData(sce)
  groups <- colData(sce)[, c("celltype3", "Case", 'Region')]
  pb <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 
  dim(pb)
  
  ## split by cluster, transform & rename columns
  pb_colData = colData(sce) %>% as.data.frame() %>%
    rownames_to_column('match') %>% 
    mutate(Pair = factor(Pair), match = paste(celltype3, Case, Region, sep = '_')) %>% 
    filter(!duplicated(match)) %>% column_to_rownames('match')
  pb_colData = pb_colData[rownames(pb),]
  
  ## make sure this PD is correct
  with(pb_colData, table(Case, celltype3, Region))
  
  ## add number of cells per aggregate
  num_cells = groups %>% as.data.frame() %>% 
    mutate(tmp = paste(celltype3, Case, Region,sep= '_')) %>% 
    pull(tmp) %>% table()
  num_cells %>% as.numeric() %>% summary()
  pb_colData$numCells = num_cells[rownames(pb_colData)]
  
  ## add the gene detection rate
  pb_colData$cdr <- scale(rowMeans(pb > 0)) 
  
  ## create SingleCellExperiment from pseudo bulk counts across all cell types and region
  (pb <- SingleCellExperiment(assays = t(pb), colData = pb_colData))
  
  ## remap case index nested inside OUD dx
  remap_case_idx = split(pb$Case, pb$DSM.IV.OUD) %>% 
    lapply(function(x){
      x = setNames(LETTERS[as.numeric(factor(x))], x)
      x[!duplicated(x)]
    }) %>% unlist()
  names(remap_case_idx) = ss(names(remap_case_idx), '\\.', 2)
  pb$CaseIdx = remap_case_idx[as.character(pb$Case)]
  
  ## check this is correct
  table(pb$DSM.IV.OUD, pb$CaseIdx)
  table(pb$celltype3, pb$CaseIdx)
  table(pb$celltype3, pb$Case)
  
  saveRDS(pb, save_pseudobulk)
} else {
  pb = readRDS(save_pseudobulk)
}

####################################################
## 2) filter pseudobulk samples that have too few cells
pb = pb[, pb$numCells > 15]
pb = pb[, pb$celltype3 != 'Mural'] # drop mural cells b/c too few
pb$celltype3 = make.names(pb$celltype3) %>% as.factor()
pb$Region = as.factor(pb$Region)
pb$Sex = as.factor(pb$Sex)
pb$numCells = as.numeric(pb$numCells)

## make an interaction term for all the combinations
pb$celltype_dx_rg_sex = interaction(pb$DSM.IV.OUD, pb$celltype3,pb$Sex, pb$Region) %>% 
  as.factor() %>% droplevels()
table(pb$celltype_dx_rg_sex)

## interaction term w/o the DX
pb$celltype_rg_sex = interaction(pb$celltype3,pb$Sex, pb$Region)  %>% 
  as.factor() %>% droplevels()
table(pb$celltype_rg_sex)

## construct design & contrast matrix regressing out the DetRate
design <- model.matrix(~ 0 + celltype_dx_rg_sex  + # term capturing the main effects
                         Age + PMI + RIN + numCells + cdr, # co-variates 
                       data = colData(pb))

## construct the null model, used in regressing out the various factors
design0 <- model.matrix(~ 0 +  celltype_rg_sex  + # term capturing the effects w/o Dx
                          Age + PMI + RIN  + numCells + cdr, # co-variates 
                        data = colData(pb))

####################################
# 3) normalization using voom-limma
y <- DGEList(counts = assays(pb)[[1]])
dim(y) # 31611   210

## filter out genes w/ low counts
A <- rowMeans(y$counts)
isexpr <- A > 5
y <- y[isexpr, , keep.lib.size = FALSE]
dim(y) # 20313   210

## filter out ribosomal genes, filter out mitochondria genes
drop.genes <- grep("^RP[SL]|^MT-",rownames(y), value = T, invert = F)
drop.genes %>% sort %>% data.frame() %>% 
  write_tsv(here(DATADIR, 'tables', 'dropped_mito_ribo_genes.tsv'))
keep.genes <- grep("^RP[SL]|^MT-",rownames(y), value = T, invert = T)
y = y[keep.genes, , keep.lib.size = FALSE]
dim(y) #  20203   210

# normalize counts
y <- calcNormFactors(y)

## voom precision weights and sample-wise quality weights normalization
v <- voomWithQualityWeights(y, design)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Case)
cor$consensus # 0.2484806

## recalculate weights after adjusting for correlated samples from same subject
v <- voomWithQualityWeights(y, design, block = colData(pb)$Case, 
          correlation = cor$consensus)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Case)
cor$consensus # 0.2465679


############################################################
# 4) Use surrogate variables to estimate unmodeled variation

## estimate the number of SVs from the adjusted
save_sva =here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sva.rds')
if(! file.exists(save_sva)){
  (n.sv = num.sv(v$E, design, method="be", seed = set.seed(1))) #21
  svobj = sva(v$E, design, design0, n.sv=n.sv, B = 20)
  saveRDS(svobj, save_sva)
} else {
  svobj = readRDS(save_sva)
}

## add the SVs to the model matrix
designSV = cbind(design, svobj$sv)
design0SV = cbind(design0, svobj$sv)

## recalculate sample quality weights after calculating the SVs
v <- voomWithQualityWeights(y, designSV, block = colData(pb)$Case, 
                            correlation = cor$consensus)
cor <- duplicateCorrelation(v, designSV, block = colData(pb)$Case)
cor$consensus # 0.08645557

save_voom = here(DATADIR, 'rdas', 'voomLimma_norm_object_N222pb.rds')
saveRDS(v, file = save_voom)

####################################################################
## 6) fit the model to get DEGs for various differences in the data
fit <- lmFit(v, designSV, block = colData(pb)$Case, correlation = cor$consensus)
fit <- eBayes(fit, robust = TRUE)

save_fit = here(DATADIR, 'rdas', 'voomLimma_diffGene_bigModelFitSVA.rds')
saveRDS(fit, file = save_fit)

save_design = here(DATADIR, 'rdas', 'bigModelFitSVA_designMatrix.rds')
saveRDS(designSV, file = save_design)

###########################################################
## 7) compute the differences b/t Dx within each cell type
celltypes=levels(factor(pb$celltype3 %>% make.names()))
designSV2 =designSV
colnames(designSV2) = make.names(colnames(designSV2))

## make the cell type contrasts
con_celltypes = sapply(setNames(celltypes, celltypes),function(cell) {
  cell = colnames(designSV2) %>% make.names() %>% str_subset(paste0('\\.', cell, '\\.'))
  OUD = cell %>% str_subset('OUD'); CTL = cell %>% str_subset('CTL')

  N_OUD = OUD %>% length(); OUD = OUD %>% paste(collapse = ' + ')
  N_CTL = CTL %>% length(); CTL = CTL %>% paste(collapse = ' + ')
  paste('(',OUD,')/',N_OUD, '-(',CTL,')/',N_CTL)
})


## proportion of each cell type
df_prop = 'data/tidy_data/tables/BU_OUD_Striatum_refined_celltype3_proportions.txt' %>% 
  read_tsv() %>% deframe()
names(df_prop) = names(df_prop) %>% make.names()
df_prop = df_prop[celltypes]

ind_neur = grepl('^D|^Int',celltypes)
ind_glia = !grepl('^D|^Int',celltypes)

## create the contrasts for OUD effect Between all cells or major classes
con_groups = c('All' = paste0('(', con_celltypes,')*', df_prop) %>% paste(collapse = ' + '), 
               'Neuron' = paste0('(', con_celltypes[ind_neur],')*', 
                                 df_prop[ind_neur]/sum(df_prop[ind_neur])) %>% paste(collapse = ' + '), 
               'Glia' =  paste0('(', con_celltypes[ind_glia],')*', 
                                df_prop[ind_glia]/sum(df_prop[ind_glia])) %>% paste(collapse = ' + '))

## refit the model based on these contrasts
cont.matrix <- makeContrasts(contrasts= c(con_groups, con_celltypes), levels=designSV2)
rownames(cont.matrix) = colnames(designSV)
fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()

## compute the DEGs from these contrasts
deg_list = lapply(setNames(colnames(cont.matrix),  names(c(con_groups, con_celltypes))), 
                  function(coef){
  topTable(coef = coef, fit =fit2, n=Inf) %>% arrange(P.Value) %>% 
  ## use SWFDR to increase power of detecting DEGs based on avg expression covariate
  ## https://pubmed.ncbi.nlm.nih.gov/30581661/
  mutate(adj.P.Val.Within =  lm_qvalue(P.Value, X=AveExpr)$q) %>%
  dplyr::select(-adj.P.Val) %>% rownames_to_column('gene') 
})

## FDR correction Between all tests
deg_list = deg_list %>% data.table::rbindlist(idcol = 'celltype') %>% 
  mutate(adj.P.Val.Between =  lm_qvalue(P.Value, X=AveExpr)$q) %>%
  split(by = 'celltype')

# FDR cutoff
sapply(deg_list, function(x) x[x$adj.P.Val.Within < 0.05,] %>% nrow())
sapply(deg_list, function(x) x[x$adj.P.Val.Between < 0.05,] %>% nrow())

# lower confidence cutoff
sapply(deg_list, function(x) x[x$P.Value < 0.01,] %>% nrow())

####################################################################
## 8) save the output of voom_limma differential state analyses
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds')
saveRDS(deg_list, save_res_fn)

tablesDir =file.path(DATADIR, 'tables'); dir.create(tablesDir, showWarnings = F)
save_res_fn2 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value)) %>% writexl::write_xlsx(save_res_fn2)

# save tables DEGs w/ P.Value < alpha up and down regulated
save_res_fn3 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.lowConfCutOff.upReg.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value) %>% filter(P.Value < 0.01, logFC > 0)) %>% 
  writexl::write_xlsx(save_res_fn3)

save_res_fn4 = here(tablesDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.lowConfCutOff.dnReg.xlsx')
deg_list %>% lapply(function(x) x %>% arrange(P.Value) %>% filter(P.Value < 0.01, logFC < 0)) %>% 
  writexl::write_xlsx(save_res_fn4)


#######################################
## 9) check out some interesting genes
sapply(deg_list, function(x){
  x %>% filter(gene %in% c('CRHR1', 'LRP8')) %>% pull(adj.P.Val.Between)
})

lapply(deg_list, function(x){
  x %>% filter(adj.P.Val.Between < 0.05) %>% pull(gene) %>% paste(collapse = ', ')
})

lapply(deg_list, function(x){
  x %>% filter(adj.P.Val.Between < 0.05) %>% pull(gene) %>% 
    str_subset('^FOX|^GAD|^CHAT$|^TH$|^KCN|^SCN|^SLC') %>% 
    paste(collapse = ', ')
})

lapply(deg_list, function(x){
  x %>% filter(adj.P.Val.Between < 0.05) %>% pull(gene) %>% 
    str_subset('^DRD|^OPR|^CADM|^CRHR') %>% paste(collapse = ', ')
})


## packages for data table processing 
library(tidyverse)
library(data.table)
library(here)
library(limma)
library(affy)
library(affycoretools)
library(annotate)
library(mouse4302.db)
library(swfdr)
library(splines)
library(org.Mm.eg.db)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/compare_to_Glud1tg_microarray'
sapply(here(DATADIR, c('rdas', 'tables')), dir.create, recursive = T, showWarnings = F)

## read in the 1-1 gene orthologs
hg38_to_mm10_genes = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/ENSEMBL_GRCh38.p13_genes_to_Orthologous_mouse_genes.tsv', show_col_types = FALSE)%>% 
  rename_with(make.names) %>% dplyr::select(c(Mouse.gene.stable.ID, Gene.name))

###################################################################
## 1) find the microarray CEL files and create the phenotype data

## first dataset from 2010 Wang et al., GSE11419
# https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-11-360
cel_fn1 = list.files(here(DATADIR,'sra'), pattern = '.CEL.gz') %>% 
  str_subset('GSM28823')

pd1 = data.frame(path = cel_fn1) %>% 
  mutate(
    Sample = basename(path) %>% ss('\\.'),
    Age = '9m', 
    Genotype = ifelse(Sample %in% paste0('GSM28823', 2:4), 'Glud1Tg', 'WT'), 
    Tissue = 'hippocampus', 
    Accession = 'GSE11419')


## second dataset from 2014 Wang et al., GSE48911
# https://bmcneurosci.biomedcentral.com/articles/10.1186/1471-2202-15-37
cel_fn2 = list.files(here(DATADIR,'sra'), pattern = '.CEL.gz') %>% 
  str_subset('GSM28823', negate = T)

pd2 = data.frame(path = cel_fn2) %>% 
  mutate(
    Sample = path %>% ss('\\_'),
    Age = path %>% ss('\\_|-|\\.', 2), 
    Genotype = ifelse(grepl('Tg', path), 'Glud1Tg', 'WT'), 
    Tissue = 'hippocampus', 
    Accession = 'GSE48911')

## combine the datasets, exclude the 
pd = bind_rows(pd1, pd2) %>% filter(Age != '10d') %>% 
  mutate(Genotype = Genotype %>% factor(c('WT', 'Glud1Tg')), 
         Age = Age %>% factor(c('4m', '9m', '14m', '20m')), 
         AgeNum = gsub('m', '', Age) %>% as.numeric(), 
         GroupName = paste(Genotype, Age)) %>% 
  column_to_rownames('Sample') 

########################################################################
## 2) import the affymatrix and create the mouse gene expression datasets
eset <- affystart(filenames = file.path( here(DATADIR, 'sra'), pd$path), 
                  phenoData = pd, groupnames = pd$GroupName, plot = F)
colnames(eset) = rownames(pd)
pData(eset) = pd

eset <- annotateEset(eset, mouse4302.db) ## annotate w/ mouse genes

## remove probes w/o ENTREZ IDs
gene_filt = fData(eset) %>% filter(!is.na(ENTREZID))
ids = mapIds(org.Mm.eg.db, keys = unique(gene_filt$ENTREZID), 
             keytype="ENTREZID", column = "ENSEMBL")
eset_filt = eset[rownames(gene_filt),]
fData(eset_filt) = cbind(fData(eset_filt), "Mouse.gene.stable.ID" = ids[gene_filt$ENTREZID])

## calculate differential expression w/ limma
mod = model.matrix(~ Accession + AgeNum + Genotype, data = pd)
arrayw <- arrayWeights(eset_filt) ## array quality weights
fit <- lmFit(eset_filt,mod, weights = arrayw) %>% eBayes()

## calculate the differential expression in the transgenic group, map to human gene
deg_tg = topTable(fit,coef="GenotypeGlud1Tg", n=Inf) %>% 
  ## use SWFDR to increase power of detecting DEGs based on avg expression covariate
  mutate(adj.P.Val =  lm_qvalue(P.Value, X=AveExpr)$q) %>% 
  filter(!duplicated(Mouse.gene.stable.ID)) %>% arrange(adj.P.Val) %>% 
  inner_join(hg38_to_mm10_genes, multiple = "all")

table(deg_tg$adj.P.Val < 0.1)

############################################################
## 3) export differential expression file from Glud1Tg mice
out_fn = here(DATADIR, 'rdas', 'Glud1Tg_microarray_reanalysis_GSE48911_GSE11419.rds')
saveRDS(deg_tg, out_fn)

xls_fn = here(DATADIR, 'tables', 'Glud1Tg_microarray_reanalysis_GSE48911_GSE11419.xlsx')
deg_tg %>% writexl::write_xlsx()

## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(DropletQC)
library(future)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

###################################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 28 cores with 500Gb
plan("multicore", workers = 28)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 150 * 1024^3)

###########################################################################
# 1) load in indvidual snRNA-seq objects and create merged Seurat projects
save_fn = list.files(here('data/raw_data/Seurat_objects'),
                     pattern = 'STARsolo_SoupX_rawCounts', full.names = T)
names(save_fn) = basename(save_fn) %>% ss('STARsolo_SoupX_rawCounts_|.rds', 2)
num_samples = length(save_fn)
objList = lapply(save_fn, readRDS)

########################################################
# 2) use Seurat reciprocal PCA to join samples together
## find integrating features
features <- SelectIntegrationFeatures(object.list = objList, nfeatures = 3000)
objList <- PrepSCTIntegration(object.list = objList, anchor.features = features)
objList <- lapply(X = objList, FUN = RunPCA, features = features, verbose = FALSE)

## select representative caudate and putamen samples for integration
## both these samples are control subjects/good number cells & depth
## and one of each sex based on:
## https://satijalab.org/seurat/articles/integration_large_datasets.html
ref = which(names(objList) %in% c('LR_RM_C1488', 'LR_RM_P1488', 
                                  'LR_RM_C13114', 'LR_RM_P1034'))

## find pair-wise anchoring cells between samples and each reference
obj.anchors <- FindIntegrationAnchors(
  object.list = objList, normalization.method = "SCT", reference = ref,
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

## merging samples together into a joint space
obj_merged <- obj.anchors %>% 
  IntegrateData(normalization.method = "SCT", dims = 1:30) %>%
  RunPCA(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>% 
  FindNeighbors(dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = 1, algorithm = 2, verbose = TRUE)


#####################################
# 3) add in patient/sample metadata
pheno = readxl::read_xlsx(here('data/tidy_data/tables/LR_RM_OUD_snRNAseq_SampleInfo.xlsx')) %>%
  rename_with(make.names) %>%
  rename_with(~ gsub("(\\.){2,}", '\\.', .x)) %>%
  rename_with(~ gsub('\\.$', '', .x)) %>%
  dplyr::rename('Infxn.Dx' = 'Lifetime.Chronic.and.or.Acute.ATOD.Infectious.or.Inflammatory.Diagnosis' ,
                'DSM.IV.SUD' = 'DSM.IV.Substance.Use.Disorder.Diagnosis.ATOD',
                'DSM.IV.Psych' = 'DSM.IV.Co.Morbid.Psychiatric.Diagnosis.ATOD', 
                'Dur.OUD' = 'Duration.of.OUD.years', 
                'Age' = 'Age.years', 
                'PMI' = 'PMI.h', 
                'pH' = 'pHa') %>%
  mutate(
    DSM.IV.OUD = ifelse(grepl('Opioid', DSM.IV.SUD), 'OUD', 'CTL') %>% factor(),
    DSM.IV.AUD = ifelse(grepl('Alcohol', DSM.IV.SUD), 'AUD', 'CTL') %>% factor(),
    DSM.IV.CUD = ifelse(grepl('Cocaine', DSM.IV.SUD), 'CUD', 'CTL') %>% factor()
  ) %>%
  column_to_rownames('Sample.ID')

## take a look at the phenotype table
head(pheno)

## look at the per-cell meta data
head(obj_merged@meta.data)

## append pt. phenotypes to single cell meta data 
obj_merged@meta.data = cbind(obj_merged@meta.data, pheno[obj_merged[[]][,'orig.ident'],])
head(obj_merged@meta.data)

###############################################
# 5) filter the lower quality cells

## Run DropletQC's identify_empty_drops function
nf.umi <- obj_merged[[]] %>% mutate(nf=dropletQC.nucFrac, umi=nCount_RNA) %>% 
  relocate(nf, umi, .before = everything())

## estimate cell, empty droplet, and damaged cell w/ DropletQC algorithms
DropletQC.ed <- nf.umi %>% 
  identify_empty_drops(nf_rescue = 0.50, umi_rescue = 1000) %>%
  relocate(cell_status, seurat_clusters, .after = umi) %>%
  dplyr::select(nf:seurat_clusters) %>%
  identify_damaged_cells(nf_sep = 0.15, umi_sep_perc = 50, verbose = FALSE)
DropletQC.ed = DropletQC.ed$df

## add Droplet QC empty droplet estimation to metadata
obj_merged$dropletQC.keep = DropletQC.ed[colnames(obj_merged), 'cell_status']

## look at which clusters should be kept by miQC per-cell fraction
obj_merged@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(miQC.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by doublet SCDS per-cell fraction
obj_merged@meta.data %>% group_by(seurat_clusters) %>%
  summarise(num = n(), prop = sum(scds.keep == 'keep') / n() ) %>% 
  arrange(prop)

## look at which clusters should be kept by both metrics
(t1 = obj_merged@meta.data %>% group_by(seurat_clusters) %>%
    summarise(num = n(), 
              numKeep = sum(scds.keep == 'keep' & miQC.keep == 'keep' & dropletQC.keep == 'cell'), 
              prop = sum(numKeep) / n() ) %>% arrange(prop))

## keep cells in the clusters that have more than 10% of OK cells
good_clusters <- t1 %>% filter(prop > 0.10) %>% pull(seurat_clusters)

## export unfiltered per-cell QC metric table
obj_merged@meta.data = obj_merged[[]] %>% relocate(dropletQC.keep, .after = 'dropletQC.nucFrac')
save_qcTale_fn = here('data/tidy_data/tables', 
                      paste0("OUD_Striatum_unfiltered_QC_table_N",num_samples,'.txt.gz'))
write_tsv(obj_merged@meta.data, save_qcTale_fn)

## subset cells to those not predicted low QC or doublet
obj_filtered = subset(obj_merged, subset = miQC.keep == "keep" & 
                        scds.keep == "keep" & dropletQC.keep == 'cell' &
                        seurat_clusters %in% good_clusters)

## recompute PCA and UMAP embedding post-filtering
obj_filtered = obj_filtered %>% RunPCA(verbose = FALSE) %>% 
  RunUMAP(reduction = "pca", dims = 1:30) %>%
  FindNeighbors(dims = 1:30, verbose = TRUE) %>% 
  FindClusters(resolution = 0.5, algorithm = 2, verbose = TRUE)


#########################################################
# 6) save projects, as h5 object for on-disk computation
save_proj_h5_fn = here('data/tidy_data/Seurat_projects', 
                       paste0("OUD_Striatum_filtered_SCT_SeuratObj_N",num_samples,'.h5Seurat'))
SaveH5Seurat(obj_filtered, filename = save_proj_h5_fn,overwrite = TRUE)

## save normalized, UMAP embedded, object for downstream analyses
save_proj_fn = here('data/tidy_data/Seurat_projects', 
                      paste0("OUD_Striatum_filtered_SCT_SeuratObj_N",num_samples,'.rds'))
saveRDS(obj_filtered, save_proj_fn)


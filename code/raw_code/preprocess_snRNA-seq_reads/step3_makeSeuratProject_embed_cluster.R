## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)
library(future)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

###################################################################
# 0) Seurat uses the future package for parallelization
## set to be parallel over 8 cores
plan("multicore", workers = 8)
options(future.rng.onMisuse = 'ignore')
options(future.globals.maxSize = 20000 * 1024^2)

###########################################################################
# 1) load in indvidual snRNA-seq objects and create merged Seurat projects
save_fn = list.files(here('data/raw_data/Seurat_objects'), 
                     pattern = 'STARsolo_SoupX_rawCounts', full.names = T)
names(save_fn) = basename(save_fn) %>% ss('STARsolo_SoupX_rawCounts_|.rds', 2)
num_samples = length(save_fn)
objList = lapply(save_fn, readRDS)

## subset cells to those not predicted low QC or doublet
objList = lapply(objList, subset, subset = miQC.keep == "keep" & scds.keep == "keep")

########################################################
# 2) use Seurat reciprocal PCA to join samples together

## find integrating features
features <- SelectIntegrationFeatures(object.list = objList, nfeatures = 3000)
objList <- PrepSCTIntegration(object.list = objList, anchor.features = features)
objList <- lapply(X = objList, FUN = RunPCA, features = features)

## find pair-wise anchoring cells for 
obj.anchors <- FindIntegrationAnchors(
  object.list = objList, normalization.method = "SCT", 
  anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)

## merging samples together into a joint space
obj_merged <- IntegrateData(anchorset = obj.anchors, normalization.method = "SCT", dims = 1:30)
obj_merged <- RunPCA(obj_merged, verbose = FALSE)
obj_merged <- RunUMAP(obj_merged, reduction = "pca", dims = 1:30)
obj_merged <- FindNeighbors(obj_merged, dims = 1:30, verbose = TRUE)
obj_merged <- FindClusters(obj_merged, resolution = 0.5, algorithm = 2, verbose = TRUE)

#####################################
# 3) add in patient/sample metadata
pheno = read.csv(here('data/tidy_data/tables/OUD1_snRNA_seq_sampleSheet.csv')) %>%
  rename_with(make.names) %>%
  rename_with(~ gsub("(\\.){2,}", '\\.', .x)) %>%
  rename_with(~ gsub('\\.$', '', .x)) %>%
  rename('Lifetime.Chronic.and.or.Acute.ATOD.Infectious.or.Inflammatory.Diagnosis' = 'Infxn.Dx',
         'DSM.IV.Substance.Use.Disorder.Diagnosis.ATOD' = 'DSM.IV.SUD',
         'DSM.IV.Co.Morbid.Psychiatric.Diagnosis.ATOD' = 'DSM.IV.Psych', 
         'Duration.of.OUD.years' = 'Dur.OUD', 
         'Age.years' = 'Age', 
         'PMI.h'= 'PMI', 
         'pHa' = 'pH') %>%
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

## export unfiltered per-cell QC metric table
save_qcTale_fn = here('data/tidy_data/tables', 
                      paste0("BU_Run1_Striatum_unfiltered_QC_table_N",num_samples,'.txt.gz'))
write_tsv(obj_merged@meta.data, save_qcTale_fn)

###################
# 4) save projects

## save as h5 object for on-disk computation
save_proj_h5_fn = here('data/tidy_data/Seurat_projects', 
                       paste0("BU_Run1_Striatum_filtered_SCT_SeuratObj_N",num_samples,'.h5Seurat'))
SaveH5Seurat(obj_merged, filename = save_proj_h5_fn,overwrite = TRUE)

## save normalized, UMAP embedded, object for downstream analyses
save_proj_fn = here('data/tidy_data/Seurat_projects', 
                      paste0("BU_Run1_Striatum_filtered_SCT_SeuratObj_N",num_samples,'.rds'))
saveRDS(obj_merged, save_proj_fn)
library(SeuratDisk)
library(Seurat)
library(tidyverse)
library(here)

## load in the filtered Seurat object
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat()

DefaultAssay(object = obj_merged) = 'RNA'
names(obj_merged[[]])

## create subset of the OUD individuals
OUD_obj = obj_merged %>% subset(DSM.IV.OUD == 'OUD')
table(OUD_obj$ID)
OUD_loom_fn = here('data/tidy_data/Seurat_projects', 
                "BU_OUD_Striatum_refined_all_SeuratObj_OUD_N12.loom")
OUD_loom <- as.loom(OUD_obj, filename = OUD_loom_fn, verbose = T)

## create a subset of the CTL individuals
CTL_obj = obj_merged %>% subset(DSM.IV.OUD == 'CTL')
table(CTL_obj$ID)
CTL_loom_fn = here('data/tidy_data/Seurat_projects', 
                   "BU_OUD_Striatum_refined_all_SeuratObj_CTL_N10.loom")
CTL_loom <- as.loom(CTL_obj, filename = CTL_loom_fn, verbose = T)




library(Seurat)
library(SeuratDisk)
library(here)

## convert the h5Seurat to h5ad file
save_merged_fn = here('data/tidy_data/Seurat_projects', 
                      "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat")


Convert(save_merged_fn, dest = "h5ad", assay = 'RNA')
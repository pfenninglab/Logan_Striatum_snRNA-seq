library(here)
library(tidyverse)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## tree specific functions
library(ggtree)
library(treeio)
library(ggtreeExtra)
library(tidytree)

## for importing python to R
library(reticulate)
pd <- import("pandas")

DATADIR='data/tidy_data/hierarchical_deg_analyses'
pickle_data <- pd$read_pickle(here(DATADIR,'tmnt_runs','LarsAIC_PCA-True',
                                   'tmnt_iteration_1.pickle'))

beta_mat = as.matrix(pickle_data$beta)
cell_assign = pickle_data$cell_assignments


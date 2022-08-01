#############################
##set working directory
setwd("/restricted/projectnb/singlecell-rl/Logan_BU_Striatum_snRNA-seq")

## Milo is available from Bioconductor (preferred stable installation)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("miloR")

## Install development version
devtools::install_github("MarioniLab/miloR", ref="devel") 

library(miloR)
library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(patchwork)

#############################

OUD_milo <- Milo(OUD_data)
OUD_milo
## class: Milo 
## dim: 29452 6875 
## metadata(0):
## assays(1): counts
## rownames(29452): ENSMUSG00000051951 ENSMUSG00000089699 ...
##   ENSMUSG00000096730 ENSMUSG00000095742
## rowData names(2): ENSEMBL SYMBOL
## colnames(6875): cell_361 cell_362 ... cell_29013 cell_29014
## colData names(17): cell barcode ... colour sizeFactor
## reducedDimNames(2): pca.corrected umap
## mainExpName: NULL
## altExpNames(0):
## nhoods dimensions(2): 1 1
## nhoodCounts dimensions(2): 1 1
## nhoodDistances dimension(1): 0
## graph names(0):
## nhoodIndex names(1): 0
## nhoodExpression dimension(2): 1 1
## nhoodReducedDim names(0):
## nhoodGraph names(0):
## nhoodAdjacency dimension(2): 1 1

OUD_milo <- buildFromAdjacency(OUD_milo, k = 30, d = 30, reduced.dim = "pca.corrected")


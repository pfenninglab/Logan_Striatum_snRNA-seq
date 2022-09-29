## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

## gene expression heatmap plots
library(Nebulosa)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'
PLOTDIR='figures/exploratory/differential_expression_analysis/plots'

## load in the DEG table
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'BU_OUD_Striatum_voom_limmaRes_bigModel_N22.rds')
res = readRDS(save_res_fn)

alpha = 0.05
df = res %>% lapply(function(x){
  x[x$p_adj < alpha, ]
}) %>% rbindlist(idcol = 'celltype')%>% 
  filter(!celltype %in% c( "DSM.IV.OUD", "RegionPutamen", "SexM" ) )

###################################
# 0) pre-set colors and cell types 
subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural/Fibroblast', 'Oligos', 'Oligos_Pre')
othertypes = c('Interneuron',  'Astrocytes', 'Endothelial', 'Microglia', 
               'Mural', 'Oligos', 'Oligos_Pre')
othertypes_col = c(carto_pal(length(othertypes) , "Vivid"))
names(othertypes_col) = othertypes

typecolors = c(subtypes_col, othertypes_col)
names(typecolors) = make.names(names(typecolors))



##################################################
# 1) load in full dataset cell type labels for plot

## read in Logan BU snRNA dataset to label transfer
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')

obj_merged$celltype3 = ifelse(grepl('Int', obj_merged$celltype3), 'Interneuron',obj_merged$celltype3)
obj_merged$celltype3 = factor(make.names(obj_merged$celltype3), names(typecolors))
table(obj_merged$celltype3)
Idents(obj_merged) = 'celltype3'



###############################################
# 2) plot some DEGs that are shared across multiple cell types
df_unique = df %>% group_by(gene) %>% 
  mutate(numCelltype = n()) %>% 
  ungroup() %>% filter(numCelltype == 1) %>% 
  filter(celltype %in% make.names(subtypes)) %>% 
  arrange(-log10(p_adj) * abs(logFC))

pdf(here(PLOTDIR, 'BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.uniqueDEGs_MSN.pdf'), width = 7.25, height = 4)
for(gene in df_unique$gene[1:50]){
  p1 = VlnPlot(obj_merged, features =gene, ncol = 1, slot = "data", pt.size = 0,
          group.by = 'celltype3', split.by = 'DSM.IV.OUD', split.plot = TRUE,
          cols = c('black', 'red'), log = F) & 
    theme(legend.position = 'none', axis.title=element_blank()) 
  print(p1)
}
dev.off()


###############################################
# 3) plot some DEGs that are cell type-specific
df_multiple = df %>% group_by(gene) %>% 
  mutate(numCelltype = n()) %>% 
  ungroup() %>% filter(numCelltype > 3) %>% 
  filter(celltype %in% make.names(subtypes)) %>% 
  arrange(-log10(p_adj) * abs(logFC) * numCelltype)

pdf(here(PLOTDIR, 'BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.multiDEGs_MSN.pdf'), width = 7.25, height = 4)
for(gene in unique(df_multiple$gene)){
  p1 = VlnPlot(obj_merged, features =gene, ncol = 1, slot = "data", pt.size = 0,
               group.by = 'celltype3', split.by = 'DSM.IV.OUD', split.plot = TRUE,
               cols = c('black', 'red'), log = F) & 
    theme(legend.position = 'none', axis.title=element_blank()) 
  print(p1)
}
dev.off()


###################################################################
# 4) plot some DEGs that are interesting for some reason or another
interesting_DEGs = c('DLGAP3', 'LRP8', 'DEPTOR', 'CRHR1')

pdf(here(PLOTDIR, 'BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.coolDEGs_MSN.pdf'), width = 7.25, height = 4)
for(gene in interesting_DEGs){
  p1 = VlnPlot(obj_merged, features =gene, ncol = 1, slot = "data", pt.size = 0,
               group.by = 'celltype3', split.by = 'DSM.IV.OUD', split.plot = TRUE,
               cols = c('black', 'red'), log = F) & 
    theme(legend.position = 'none', axis.title=element_blank()) 
  print(p1)
}
dev.off()


df %>% filter(gene == 'CRHR1') %>% as.data.frame()
df %>% filter(gene == 'LRP8') %>% as.data.frame()
df %>% filter(gene == 'CELF4') %>% as.data.frame()
df %>% filter(gene == 'SH3BP4') %>% as.data.frame()

df %>% filter(grepl('HIF', gene)) %>% as.data.frame()
df %>% filter(grepl('RELN', gene)) %>% as.data.frame()
df %>% filter(grepl('ADAMTSL2', gene)) %>% as.data.frame()
df %>% filter(grepl('OSER1', gene)) %>% as.data.frame()


###############################################
# 5) plot some DEGs that are shared across multiple cell types
df_unique = df %>% group_by(gene) %>% 
  mutate(numCelltype = n()) %>% 
  ungroup() %>% filter(numCelltype == 1) %>% 
  filter(celltype %in% make.names(othertypes)) %>% 
  arrange(log10(p_adj) * abs(logFC) * AveExpr)

pdf(here(PLOTDIR, 'BU_OUD_Striatum_voomLimmaRes_N22.bycelltype3.uniqueDEGs_glia.pdf'), width = 7.25, height = 4)
for(gene in df_unique$gene[1:50]){
  p1 = VlnPlot(obj_merged, features =gene, ncol = 1, slot = "data", pt.size = 0,
               group.by = 'celltype3', split.by = 'DSM.IV.OUD', split.plot = TRUE,
               cols = c('black', 'red'), log = F) & 
    theme(legend.position = 'none', axis.title=element_blank()) 
  print(p1)
}
dev.off()




#########################################################
## plot the main UMAP of the refined cell type labels ###
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_UMAP_all.celltype3byDx.pdf'), width = 10.5, height = 8)
DimPlot(object = obj_merged, reduction = "umap", label.size = 3,  
        pt.size = 2, group.by = 'celltype3',  split.by = c('DSM.IV.OUD'),
        cols = typecolors) +
  guides(color = guide_legend(nrow = 4, override.aes= list(size = 1)))+ 
  theme(legend.position = 'bottom', legend.spacing.x = unit(1.0, 'cm'))
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_UMAP_all.celltype3byRegion.pdf'), width = 10.5, height = 8)
DimPlot(object = obj_merged, reduction = "umap", label.size = 3,  
        pt.size = 2, group.by = 'celltype3',  split.by = c('Region'),
        cols = typecolors) +
  guides(color = guide_legend(nrow = 4, override.aes= list(size = 1)))+ 
  theme(legend.position = 'bottom', legend.spacing.x = unit(1.0, 'cm'))
dev.off()





## Neuron vs. Glia markers
markerGenes1 <- c('RBFOX3', 'LHX6', 'AQP4', 'CX3CR1', 'PDGFRA', 'MOG' )
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_UMAP_all.allMarkers.pdf'), width = 7.25, height = 4)
p1 = plot_density(obj_merged, features = markerGenes1,  reduction = "umap") +
  plot_layout(nrow = 2, guides = "auto") &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 10))
p1
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_all.allMarkers.pdf'), width = 7.25, height = 5)
VlnPlot(obj_merged, features =c(markerGenes1), slot = "data",  pt.size = 0,
        cols = c('MSNs' = 'gray', othertypes_col)) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
dev.off()


## MSN subtypes
markMSN1 = c('Drd1','Tac1','Reln') %>% toupper()# D1 markers
markMSN2 = c('Drd2','Adora2a','Penk')%>% toupper() # D2 markers
markMSN3 = c('Foxp2', 'Rxfp1', 'Casz1')%>% toupper() # D1/2H markers

pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_UMAP_all.MSNclusterMarkers.pdf'), width = 7.25, height = 3)
p2 = plot_density(obj_merged, slot = 'data', features = markMSN1,  reduction = "umap", joint= T) &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 8))
p3 = plot_density(obj_merged, slot = 'data', features = markMSN2,  reduction = "umap", joint= T) &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 8))
p4 = plot_density(obj_merged, slot = 'data', features = markMSN3,  reduction = "umap", joint= T) &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 8))
p2 + plot_layout(nrow = 1) & theme(legend.position = 'bottom')
p3 + plot_layout(nrow = 1) & theme(legend.position = 'bottom')
p4 + plot_layout(nrow = 1) & theme(legend.position = 'bottom')
dev.off()


## MSN compartments
markMSN4 = c('STXBP6', 'SEMA3E', 'EPHA4', 'GDA') # Matrix markers
markMSN5 =c( 'PDYN', 'OPRM1', 'KHDRBS3', 'KCNIP1') # Striosome markers

pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_UMAP_all.MSNcompMarkers.pdf'), width = 7.25, height = 3)
p5 = plot_density(obj_merged, slot = 'data', features = markMSN4,  reduction = "umap") 
p6 = plot_density(obj_merged, slot = 'data', features = markMSN5,  reduction = "umap", pal = 'inferno')
p5 + plot_layout(nrow = 1) & theme_classic(base_size = 7) & 
  theme(plot.title = element_text(size = 8)) & theme(legend.position = 'bottom') 
p6 + plot_layout(nrow = 1) & theme_classic(base_size = 7) & 
  theme(plot.title = element_text(size = 8)) & theme(legend.position = 'bottom') 
dev.off()




## Striatal interneuron markers, PMID: 30134177, Munoz-Manchado et al.
markINT1 =c( 'LHX6', 'SST', 'NPY', 'CHAT', 'PTHLH','PVALB', 'TH','TRH') 
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_UMAP_all.InterneuronMarkers.pdf'), width = 7.25, height = 3)
p7 = plot_density(obj_merged, slot = 'data', features = markINT1,  reduction = "umap") 
p7 + plot_layout(nrow = 2) & theme_classic(base_size = 7) & 
  theme(plot.title = element_text(size = 8))

# p8 = plot_density(obj_merged, slot = 'data', features = markINT2,  reduction = "umap") 
# p8 + plot_layout(nrow = 2) & theme_classic(base_size = 7) & 
#   theme(plot.title = element_text(size = 8))
dev.off()







## plot the opioid receptors and ligands
markOPR1 = c('OPRD1', 'OPRM1', 'OPRK1') # opioid receptor genes
markOPR2 = c( 'PENK', 'PDYN') # endorphin peptide genes
obj_merged$celltype1 = factor(obj_merged$celltype1, levels = c('MSNs', othertypes)) 
Idents(obj_merged) = 'celltype1'

pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_all.opioids.pdf'), width = 7.25, height = 5)
VlnPlot(obj_merged, features =c(markOPR1, markOPR2), ncol = 3, slot = "data", pt.size = 0,
        fill.by = 'celltype1', group.by = 'celltype1', 
        cols = c('MSNs' = 'gray', othertypes_col), log = T) & 
  theme(legend.position = 'none', axis.title=element_blank()) & 
  coord_flip() 
dev.off()

pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_all.opioidsByDxSUD.pdf'), width = 7.25, height = 5)
VlnPlot(obj_merged, features =c(markOPR1, markOPR2), slot = "data", pt.size = 0, 
        ncol = 3, group.by = 'celltype1', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), 
        split.plot = TRUE, log = T) & 
  theme(axis.title=element_blank()) & coord_flip() 
dev.off()


pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_all.opioidsByDxSUD.pdf'), width = 7.25, height = 5)
VlnPlot(obj_merged, features =c(markOPR1, markOPR2), slot = "data", pt.size = 0, 
        ncol = 3, group.by = 'celltype1', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), 
        split.plot = TRUE, log = T) & 
  theme(axis.title=element_blank()) & coord_flip() 
dev.off()




#####################################
# 2) plot genes at MSN subtype level
## read in Logan BU snRNA dataset to label transfer
obj_msn = here('data/tidy_data/Seurat_projects', 
               "BU_OUD_Striatum_subsetMSN_SCT_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA')
obj_msn$celltype3 = factor(obj_msn$celltype3 , c(subtypes, othertypes))
Idents(obj_msn) = 'celltype3'

## plot by direct and indirect MSN markers, 
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_msn.MSNclusterMarkers.pdf'), width = 7.25, height = 2.5)
VlnPlot(obj_msn, features =c(markMSN1), slot = "data", ncol = 3, pt.size = 0, 
        group.by = 'celltype3', cols = typecolors) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
VlnPlot(obj_msn, features =c(markMSN2), slot = "data", ncol = 3, pt.size = 0, 
        group.by = 'celltype3', cols = typecolors) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
VlnPlot(obj_msn, features =c(markMSN3), slot = "data", ncol = 3, pt.size = 0,
        group.by = 'celltype3', cols = typecolors) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
dev.off()


## plot by direct and indirect MSN markers, stratify by Dx
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_msn.MSNclusterMarkersByDxSUD.pdf'), width = 7.25, height = 2.5)
VlnPlot(obj_msn, features =c(markMSN1), slot = "data", ncol = 3, pt.size = 0, 
        group.by = 'celltype3', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), split.plot = TRUE) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
VlnPlot(obj_msn, features =c(markMSN2), slot = "data", ncol = 3, pt.size = 0, 
        group.by = 'celltype3', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), split.plot = TRUE) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
VlnPlot(obj_msn, features =c(markMSN3), slot = "data", ncol = 3, pt.size = 0,
        group.by = 'celltype3', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), split.plot = TRUE) & 
  theme(legend.position = 'none', axis.title=element_blank()) 
dev.off()


## plot by patch/matrix markers
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_msn.MSNcompMarkers.pdf'), width = 7.25, height = 5)
VlnPlot(obj_msn, features =c(markMSN4, markMSN5), slot = "data", ncol = 4, pt.size = 0, 
        group.by = 'celltype3', cols = typecolors) & 
  theme(legend.position = 'none', axis.title=element_blank())
dev.off()


pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_msn.MSNcompMarkersByDxSUD.pdf'), width = 7.25, height = 5)
VlnPlot(obj_msn, features =c(markMSN4, markMSN5), slot = "data", ncol = 4, pt.size = 0, 
        group.by = 'celltype3', split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), split.plot = TRUE) & 
  theme(legend.position = 'none', axis.title=element_blank())
dev.off()


## plot by opioid receptors or endogenous peptides
pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_msn.opioids.pdf'), width = 7.25, height = 5)
VlnPlot(obj_msn, features =c(markOPR1, markOPR2), slot = "data", cols = subtypes_col, pt.size = 0) & 
  theme(legend.position = 'none', axis.title=element_blank()) & coord_flip() 
dev.off()


pdf(here(PLOTDIR, 'plots', 'BU_OUD_Striatum_Viol_msn.opioidsByDxSUD.pdf'), width = 7.25, height = 5)
VlnPlot(obj_msn, features =c(markOPR1,markOPR2), slot = "data", ncol = 3, pt.size = 0,
        split.by = 'DSM.IV.OUD', cols = c('gray', 'red'), split.plot = TRUE) & 
  theme(legend.position = 'none', axis.title=element_blank()) & coord_flip() 
dev.off()



## conda activate r4
## packages for data table processing 
library(here)
library(tidyverse)
library(RColorBrewer)
library(rcartocolor)
library(ggpubr)

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

PLOTDIR='figures/exploratory/preprocess_snRNA-seq_reads'

##################################################
# 1) load in full dataset cell type labels for plot

## load the unfiltered QC table
qc_df = here('data/tidy_data/tables',
             paste0("BU_Run1_Striatum_unfiltered_QC_table_N4.txt.gz")) %>%
  read_tsv(show_col_types = FALSE) %>%
  mutate(toKeep = case_when(
    scds.keep == 'doublet' ~ 'doublet', 
    dropletQC.keep == 'empty_droplet' ~ 'empty droplet',
    miQC.keep == 'discard' ~ 'high mito rate', 
    TRUE ~ 'keep'), 
    toKeep = factor(toKeep, c('doublet', 'empty droplet', 'high mito rate', 'keep')))

## plot number of estimated miQC compromised cells
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_QC_all.perSampleQC.pdf'), 
    width = 7.25, height = 4, onefile = F)
p1 = qc_df %>% ggplot(aes(x = orig.ident, fill = toKeep)) + 
  geom_bar(stat = 'count')  + ylab('Number of cells') + 
  scale_fill_carto_d(palette = 'ArmyRose') + 
  theme_classic() + xlab('Sample')
p2 =  qc_df %>% ggplot(aes(x = orig.ident, fill = toKeep)) + 
  geom_bar(stat = 'count', position = 'fill') + 
  scale_fill_carto_d(palette = 'ArmyRose') + 
  theme_classic() + ylab('Proportion of cells') + 
  xlab('Sample')

ggarrange(p1, p2, labels = c("A", "B"), 
          common.legend = TRUE, nrow = 1, legend="bottom")
dev.off()

## load in the filtered Seurat object
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_Run1_Striatum_filtered_SCT_SeuratObj_N4.h5Seurat") %>% LoadH5Seurat() 

## plot number of estimated miQC compromised cells
qc_df2 = obj_merged[[]]
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_QC_all.qcViolinPlots.pdf'), 
    width = 7.25, height = 4, onefile = F)
p1 = ggviolin(qc_df2, x = "orig.ident", y = "percent.mt", add = "boxplot",
              palette = c("gray", "red"), fill = 'DSM.IV.OUD') + 
  ylab('Percent Mitochondrial reads') +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1), 
        axis.title.x=element_blank())
p2 = ggviolin(qc_df2, x = "orig.ident", y = "dropletQC.nucFrac", add = "boxplot",
              palette = c("gray", "red"), fill = 'DSM.IV.OUD') +
  ylab('Percent nuclear reads') +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1), 
        axis.title.x=element_blank())
  
ggarrange(p1, p2, labels = c("A", "B"), 
          common.legend = TRUE, nrow = 1, legend="bottom")
dev.off()


## plot the per-cell QC metric data on the UMAP
obj_merged$nLogUMI = log10(obj_merged$nCount_RNA)
DefaultAssay(obj_merged) = 'RNA'

## Per Cell type QC values
QC_features = c('percent.mt',  'scds.hybrid_score', 'dropletQC.nucFrac', 'nLogUMI')
pdf(here(PLOTDIR, 'plots', 'BU_Run1_Striatum_QC_all.perCellQC.pdf'), width = 7.25, height = 4)
FeaturePlot(object = obj_merged, reduction = "umap", feature = QC_features)
plot_density(obj_merged, features = QC_features,  reduction = "umap") +
  plot_layout(nrow = 2, guides = "auto") &
  theme_classic(base_size = 7) & theme(plot.title = element_text(size = 10))
dev.off()










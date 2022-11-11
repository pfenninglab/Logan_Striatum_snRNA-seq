ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(ggsci)

## main Seurat package snRNA-seq pacakges
library(Seurat)
library(SeuratDisk)
library(future)

library(here)

DATADIR='data/raw_data'

## make for this subdirs
PLOTDIR='figures/explanatory/figure1_the_experiment'
here(PLOTDIR, c('plots', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F)

###################################
# 0) pre-set plotting aesthetics
in2mm<-25.4
my_theme = theme_classic(base_size = 6)

###################################
# 1) read in Logan snRNA dataset to plot
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA') 
names(obj_merged[[]] )
table(obj_merged$celltype3 )

###################################
# 2) plot the sample-wise QC metrics
df1 = obj_merged[[]] %>% 
  group_by(ID, Region, DSM.IV.OUD) %>% 
  summarise(`# nuclei` = n())

## load the unfiltered QC table
df2 = here('data/tidy_data/tables',
             paste0("BU_OUD_Striatum_unfiltered_QC_table_N22.txt.gz")) %>%
  read_tsv(show_col_types = FALSE) %>%
  group_by(ID, Region, DSM.IV.OUD) %>% 
  summarise(numTotal = n()) %>% 
  inner_join(df1) %>% 
  mutate(`% passing QC` = `# nuclei`/numTotal * 100 )%>% 
  pivot_longer(cols = c(`# nuclei`, `% passing QC`, numTotal), 
               names_to = 'metric')

fig1_violin_propCells_per_sample_fn = 
  here(PLOTDIR, 'plots', 'fig1_violin_Ncell_propCells_per_sample.pdf')

pdf(fig1_violin_propCells_per_sample_fn, width = 40/in2mm, height =  50/in2mm)
ggplot(df2, aes(x = DSM.IV.OUD, y = value, fill = DSM.IV.OUD)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 21, size = 1) + 
  facet_grid(metric ~ Region, scales = 'free_y') + 
  my_theme + scale_fill_npg() + 
  scale_y_continuous(limits =c(0, NA), labels = scales::comma) + 
  ylab('per-sample QC metric') +
  theme(legend.position = 'none', axis.title.x = element_blank())
dev.off()



###################################
# 3) plot the avg. cell-wise QC metrics
df3 = obj_merged[[]] %>% 
  group_by(ID, Region, DSM.IV.OUD) %>% 
  summarise(`Avg. Genes` = mean(nFeature_RNA), 
            `Avg. UMI` = mean(nCount_RNA)) %>% 
  pivot_longer(cols = c(`Avg. Genes`, `Avg. UMI`), 
               names_to = 'metric')

fig1_violin_avgUMI_per_sample_fn = 
  here(PLOTDIR, 'plots', 'fig1_violin_avgUMI_avgGene_per_sample.pdf')

pdf(fig1_violin_avgUMI_per_sample_fn, width = 40/in2mm, height =  50/in2mm)
ggplot(df3, aes(x = DSM.IV.OUD, y = value, fill = DSM.IV.OUD)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(pch = 21, size = 1) + 
  facet_grid(metric ~ Region, scales = 'free_y') + 
  my_theme + scale_fill_npg() +
  scale_y_continuous(trans='log10', labels = scales::comma) + 
  ylab('per-nuclei QC metric') +
  theme(legend.position = 'none', axis.title.x = element_blank())
dev.off()


############################
# 4) export the average su
stable1_per_sample_QC_fn = 
  here(PLOTDIR, 'tables', 's1.3_table_post_sequencing_qc_per_sample.xlsx')
df4 = bind_rows(df2, df3) %>% 
  pivot_wider(names_from = 'metric', values_from = 'value') %>% 
  writexl::write_xlsx(stable1_per_sample_QC_fn)


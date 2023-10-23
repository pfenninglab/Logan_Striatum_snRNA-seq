ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(tidymodels)
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

#######################################################
# 0) Seurat uses the future package for parallelization
plan("sequential")
options(future.globals.maxSize = 100 * 1024^3)
options(future.rng.onMisuse = 'ignore')

in2mm<-25.4; alpha = 0.05
my_theme = theme_classic(base_size = 5)

subtypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome','D1/D2-Hybrid')
subtypes_col = c('#1f78b4', '#a6cee3', '#e31a1c', '#fb9a99',  '#6a3d9a')
names(subtypes_col) = subtypes %>% make.names()

othertypes = c('Int-CCK', 'Int-PTHLH','Int-SST', 'Int-TH', 
               'All', 'Neuron', 'Glia',
               'Astrocytes', 'Endothelial', 'Microglia', 'Mural', 'Oligos', 
               'Oligos_Pre', 'Interneurons')
othertypes_col = c(carto_pal(4, "Safe"), 
                   mypal = pal_npg('nrc')(3),
                   carto_pal(length(othertypes) -7 , "Vivid"))
names(othertypes_col) = othertypes
othertypes_col = othertypes_col[-c(1:4)]
typecolors = c(subtypes_col, othertypes_col[c(10, 4:9)])


###################################
# 1) read in Logan snRNA dataset to plot
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat(assay = 'RNA') 
names(obj_merged[[]] )
table(obj_merged$celltype3 )

obj_merged$celltype4 = ifelse(grepl('Int', obj_merged$celltype3), 
                              'Interneurons', obj_merged$celltype3) %>% make.names()
table(obj_merged$celltype4 )

###################################
# 2) plot the sample-wise QC metrics
df1 = obj_merged[[]] %>% 
  group_by(ID, Region, DSM.IV.OUD) %>% 
  summarise(`# nuclei` = n())

## load the unfiltered QC table
df2 = here('data/tidy_data/tables',
             paste0("BU_OUD_Striatum_unfiltered_QC_table_N22.txt.gz")) %>%
  read_tsv(show_col_types = FALSE) %>%
  group_by(ID, Region, DSM.IV.OUD, Sex) %>% 
  summarise(numTotal = n()) %>% 
  inner_join(df1) %>% 
  mutate(`% passing QC` = `# nuclei`/numTotal * 100 )%>% 
  pivot_longer(cols = c(`# nuclei`, `% passing QC`, numTotal), 
               names_to = 'metric')

## get descriptive numbers
df2 %>% group_by(Region, metric) %>% 
  summarise(num = mean(value),
            se = sd(value)/sqrt(n()))

## t-test of the differences
df2 %>% nest(data = -c(Region, metric)) %>% 
  mutate(
    test = map(data, ~ t.test(value~DSM.IV.OUD, data = .x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  select(-data, -test) %>% as.data.frame() %>% 
  writexl::write_xlsx(  
    here(PLOTDIR, 'tables', 'fig1_violin_Ncell_propCells_per_sample.xlsx')
)


## plot the 
fig1_violin_propCells_per_sample_fn = 
  here(PLOTDIR, 'plots', 'fig1_violin_Ncell_propCells_per_sample.pdf')

pdf(fig1_violin_propCells_per_sample_fn, width = 40/in2mm, height =  50/in2mm)
ggplot(df2, aes(x = DSM.IV.OUD, y = value, fill = DSM.IV.OUD)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = Sex), size = 1) + 
  facet_grid(metric ~ Region, scales = 'free_y') + 
  my_theme + scale_fill_npg() + 
  scale_y_continuous(limits =c(0, NA), labels = scales::comma) + 
  ylab('per-sample QC metric') +
  theme(legend.position = 'none', axis.title.x = element_blank())
dev.off()


###################################
# 3) plot the avg. cell-wise QC metrics
df3 = obj_merged[[]] %>% 
  group_by(ID, Region, DSM.IV.OUD, Sex) %>% 
  summarise(`Avg. Genes` = mean(nFeature_RNA), 
            `Avg. UMI` = mean(nCount_RNA)) %>% 
  pivot_longer(cols = c(`Avg. Genes`, `Avg. UMI`), 
               names_to = 'metric')

## get descriptive numbers
df3 %>% group_by(metric) %>% 
  summarise(num = mean(value),
            se = sd(value)/sqrt(n()))

## t-test of the differences
df3 %>% nest(data = -c(Region, metric)) %>% 
  mutate(
    test = map(data, ~ t.test(value~DSM.IV.OUD, data = .x)), # S3 list-col
    tidied = map(test, tidy)
  ) %>% 
  unnest(cols = tidied) %>% 
  select(-data, -test) %>% as.data.frame() %>% 
  writexl::write_xlsx(  
    here(PLOTDIR, 'tables', 'fig1_violin_avgUMI_avgGene_per_sample.xlsx')
  )

fig1_violin_avgUMI_per_sample_fn = 
  here(PLOTDIR, 'plots', 'fig1_violin_avgUMI_avgGene_per_sample.pdf')

pdf(fig1_violin_avgUMI_per_sample_fn, width = 40/in2mm, height =  50/in2mm)
ggplot(df3, aes(x = DSM.IV.OUD, y = value, fill = DSM.IV.OUD)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = Sex), size = 1) + 
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



###################################
# 5) plot the avg. cell-wise QC metrics
df4 = obj_merged[[]] %>% 
  group_by(Case, DSM.IV.OUD, Sex) %>% 
  summarise(`Avg. Genes` = mean(nFeature_RNA), 
            `Avg. UMI` = mean(nCount_RNA), 
            `num Cells` = n()) %>% 
  pivot_longer(cols = c(`Avg. Genes`, `Avg. UMI`, `num Cells`), 
               names_to = 'metric')

## get descriptive numbers
df4 %>% group_by(metric) %>% 
  summarise(num = mean(value),
            se = sd(value)/sqrt(n()))



#############################################
# 6) plot the avg. cell type wise QC metrics
df5 = obj_merged[[]] %>% distinct() %>% 
  group_by(Case, DSM.IV.OUD, Region, Sex, celltype4) %>% 
  summarise(`Avg. Genes` = mean(nFeature_RNA), 
            `Avg. UMI` = mean(nCount_RNA)) %>% 
  pivot_longer(cols = c(`Avg. Genes`, `Avg. UMI`), 
               names_to = 'metric') %>% 
  mutate(
    celltype4 = factor(celltype4, levels = names(typecolors)),
    cell_class = case_when(
    celltype4 %in% c(names(subtypes_col), 'Interneurons') ~ 'Neuron',
    TRUE ~ 'Glia'
  ), cell_class = factor(cell_class, levels = c('Neuron', 'Glia')))


df5 %>% group_by(metric, cell_class) %>% 
  summarise(mean = mean(value),
            se = sd(value)/sqrt(n())) 

fig1_violin_avgUMI_per_sample_fn = 
  here(PLOTDIR, 'plots', 'sfig1_violin_avgUMI_avgGene_per_sample_per_celltype.pdf')

pdf(fig1_violin_avgUMI_per_sample_fn, width = 120/in2mm, height =  80/in2mm)
ggplot(df5, aes(x = celltype4, y = value, fill = celltype4)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(aes(shape = Sex), size = 1) + 
  facet_grid(metric ~ cell_class, scales = 'free', space = 'free_x') + 
  my_theme + 
  scale_fill_manual(values = typecolors) +
  scale_shape_manual(values = c(21, 22)) +
  scale_y_continuous(labels = scales::comma, limits = c(0,NA)) + 
  ylab('Average per-nuclei QC metric across cell types') +
  theme(legend.position = 'none', axis.title.x = element_blank())
dev.off()


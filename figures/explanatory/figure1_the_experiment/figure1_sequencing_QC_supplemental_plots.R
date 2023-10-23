ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(tidyverse)
library(rcartocolor)
library(ggsci)
library(cowplot)
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
# 1) grab the subject phenotype data
pheno =readRDS(here('data/tidy_data/tables/LR_RM_OUD_snRNAseq_SampleInfo.rds')) %>% 
  rownames_to_column('Sample') %>% 
  filter(!ID %in% c('C-13291', 'C-612'))


###################################
# 2) load the mapping to genome QC table
df1 = here('figures/exploratory/explore_qc_data/tables',
           'Logan-OUD-snRNAseq.mapping_QC_stats.BU_Run0-7_20220517.xlsx') %>%
  readxl::read_xlsx() 
df1 = df1[, c(T, T, apply(df1[,-c(1:2)], 2, function(x) var(x) > 0.05))]
df1 = inner_join(x = pheno, y = df1) %>% 
  pivot_longer(cols = -all_of(names(pheno)), names_to = 'metric')

p1 = ggplot(df1, aes(x = Region, y = value, fill = DSM.IV.OUD)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(pch = 21, size = 1, position = position_jitterdodge()) + 
  my_theme + scale_fill_npg() + 
  facet_wrap(~metric , scales = 'free_y', labeller = label_wrap_gen(width = 28), ncol = 5 ) + 
  scale_y_continuous(limits =c(0, NA), labels = scales::comma) + 
  ylab('per-sample mapping QC metric') +
  theme(legend.position = 'bottom', axis.title.x = element_blank())


###################################
# 3) grab the avg. cell-wise QC metrics
df2 = here('figures/exploratory/explore_qc_data/tables',
           'Logan-OUD-snRNAseq.summary_QC_stats.BU_Run0-7_20220517.xlsx') %>%
  readxl::read_xlsx() %>% dplyr::select(-Notes)
df2 = df2[, c(T, T, apply(df2[,-c(1:2)], 2, function(x) var(x) > 0.05))]
df2 = inner_join(x = pheno, y = df2) %>% 
  pivot_longer(cols = -all_of(names(pheno)), names_to = 'metric')

p2 = ggplot(df2, aes(x = Region, y = value, fill = DSM.IV.OUD)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_point(pch = 21, size = 1, position = position_jitterdodge()) + 
  my_theme + scale_fill_npg() + 
  facet_wrap(~metric , scales = 'free_y', labeller = label_wrap_gen(width = 28), ncol = 5) + 
  scale_y_continuous(limits =c(0, NA), labels = scales::comma) + 
  ylab('per-sample single-cell QC metric') +
  theme(legend.position = 'bottom', axis.title.x = element_blank())


###################################
# 4) make the combined QC plot
fig1_alignment_qc_fn = 
  here(PLOTDIR, 'plots', 's1.4_figure_alignment_qc.pdf')

pdf(fig1_alignment_qc_fn, width = 180/in2mm, height =  210/in2mm)
plot_grid( p1, p2 + theme(legend.position="none"),
  align = 'vh', labels = c("A", "B"), hjust = -1, ncol = 1, 
  rel_heights = c(2, 1.35))
dev.off()

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)

## packages for data table processing/plotting
library(data.table)
library(tidyverse)
library(tidymodels)
library(here)
library(ggh4x)
library(cowplot)
library(ggpubr)

DATADIR='/projects/pfenninggroup/singleCell/McLean_chronic_opioid_monkey_snRNA-seq/data/tidy_data/differential_expression'

in2mm<-25.4
my_theme = theme_classic(base_size = 6)

## make for this subdirs
PLOTDIR='figures/exploratory/compare_McLean_monkey_snRNA-seq'
here(PLOTDIR, c('plots/RRHO2', 'tables', 'rdas')) %>% sapply(dir.create, showWarnings = F, recursive = T)

##################################################################
# 1) grab the OUD DEGs split by the DEGs to plot in neurons or glia
alpha = 0.05
human_deg = here('data/tidy_data/differential_expression_analysis', 'rdas', 'OUD_Striatum_voom_limma_bigModelSVA_N22.celltype.rds') %>% 
  readRDS() %>% bind_rows() %>% rename('celltype_hs' = 'celltype') %>% 
  rename_if(is.numeric, ~ paste0(., '_hs')) %>% 
  filter(P.Value_hs < alpha)

###################################################################
# 2) load in the DEGs from monkey snRNA-seq NAc exposed to opioids
all_degs_fn = list.files(here(DATADIR, 'rdas'), pattern = 'celltype.rds', full.names = T, recursive = F) 

monkey_deg_versions = data.frame(path = all_degs_fn) %>% 
  mutate(deg_version = dirname(path) %>% dirname() %>% basename(), 
         countsFilter = ss(deg_version, '_', 2), covariate = 'Pair', 
         sva = basename(path) %>% ss('_', 3), 
         theCall = paste(countsFilter, covariate, sva, sep = '_')) %>% 
  filter(deg_version != 'Result_batch0_with_mistakes')

monkey_deg_paths = monkey_deg_versions %>% dplyr::select(c(theCall, path)) %>% deframe()
monkey_deg_name = monkey_deg_versions %>% pull(theCall) %>% set_names()

drop_combo = c("D1.ICj : D1.D2.Hybrid", "D1.ICj : D1.Matrix",
               "D1.ICj : D1.Striosome", "D1.NUDAP : D1.Matrix",
               "D1.NUDAP : D1.Striosome", "D1.Shell.OT : D1.D2.Hybrid",
               "Oligos_Pre : Oligos", "Oligos : Oligos_Pre")

#############################################################################
# 2) read in the monkey DEGs from each version of the deg computation
deg_name = "expression_Pair_sva" 
monkey_deg = readRDS(monkey_deg_paths[deg_name])  %>% 
bind_rows() %>% rename_if(is.numeric, ~ paste0(., '_rm')) %>% 
rename('celltype_rm' = 'celltype')

full_deg = inner_join(human_deg, monkey_deg, multiple = "all") %>% 
  mutate(combo = paste(celltype_rm, celltype_hs, sep = ' : ')) %>% 
  filter(celltype_rm %>% str_sub(1,3) == celltype_hs %>% str_sub(1,3)) %>% 
  filter(!combo %in% c(drop_combo))  %>% 
  mutate(replicated = case_when(
    P.Value_rm < alpha & adj.P.Val.Between_hs < alpha ~ T, 
    P.Value_rm > alpha & adj.P.Val.Between_hs < alpha ~ F, 
  ), 
  replicated = factor(replicated)) %>% 
  filter(!is.na(replicated)) %>% 
  mutate(group = case_when(grepl('All|Glia|Neuron', combo) ~ 'Class', 
                           grepl('^D|^In', combo) ~ 'Neuron subtype', 
                           T ~ 'Glia type'))
table(full_deg$replicated)

################################################
# 3) plot the number of replicated human DEGs 
num_replicated_df = full_deg %>% 
  group_by(celltype_rm, celltype_hs, group, combo) %>% count(replicated) %>% 
  arrange(combo, rev(replicated)) 
## make the plot, number and percent
plot_fn = here(PLOTDIR,'plots', paste0('figureNN_compare_McLean_monkey_snRNA-seq.percentReplication.', deg_name, '.pdf' ))
pdf(plot_fn, width = 180/in2mm, height = 80/in2mm)

## number of replicated DEGs
p1 = ggplot(num_replicated_df, aes(x = combo, y = n, fill = replicated)) + 
  geom_bar(stat="identity") + coord_flip() + 
  geom_text(aes(label = n), size = 2, hjust = 'inside') +
  facet_grid(group ~ ., scales = 'free_y', space = 'free') +
  scale_x_discrete(limits = rev) + my_theme + 
  scale_fill_brewer(palette = 'Paired', name = 'Morphine\nvs. Control\nP < 0.05') +
  theme(legend.position = 'none') + 
  ylab('Number of human DEGs expressed in monkey cell type') + 
  xlab('Cell type')

## calculate and plot the percent replicated DEGs
p2 = num_replicated_df %>% group_by(combo) %>% 
  mutate(pct= n/sum(n) * 100) %>% 
  ggplot(aes(x = combo, y = pct, fill = replicated)) + 
  geom_bar(stat="identity") + coord_flip() + 
  geom_text(aes(label = signif(pct, 2)), size = 2, hjust = 'inside') +
  facet_grid(group ~ ., scales = 'free_y', space = 'free') +
  scale_x_discrete(limits = rev) + my_theme + 
  scale_fill_brewer(palette = 'Paired', name = 'Morphine\nvs. Control\nP < 0.05') +
  theme(legend.position = 'right', axis.text.y = element_blank()) + 
  ylab('Percent of human DEGs expressed in monkey cell type') + 
  xlab('Cell type')

plot_grid(p1, p2, labels = c('A', 'B'), nrow =1, label_size = 6, rel_widths = c(1, 1))
dev.off()

################################################
# 4) plot the number of replicated human DEGs 
cor_lfc = full_deg %>% 
  nest(data = -c(combo, group)) %>% 
  mutate(
    n = map_int(data, nrow),
    test = map(data, ~ cor.test(.x$logFC_hs, .x$logFC_rm)),
    tidied = map(test, tidy)
  ) %>% unnest(cols = tidied) %>% 
  select(-data, -test) %>% arrange(p.value)

# plot the RRHO2 grouped all together with one scale
plot_fn = here(PLOTDIR,'plots', paste0('figureNN_compare_McLean_monkey_snRNA-seq.lfcCorrelation.', deg_name, '.pdf' ))
pdf(plot_fn, width = 180/in2mm, height = 240/in2mm)

ggplot(full_deg, aes(x = logFC_hs, y = logFC_rm)) + 
  geom_hline(yintercept = 0, color = 'black') + 
  geom_vline(xintercept = 0, color = 'black') + 
  geom_point(pch = 21, aes(fill = replicated, alpha = replicated)) + 
  geom_smooth(method = 'lm') + 
  stat_cor(size = 2) +
  facet_nested_wrap( ~ group + combo, ncol = 3) + 
  scale_fill_brewer(palette = 'Paired', name = 'Morphine vs. Control P < 0.05') +
  scale_alpha_manual(values = c(0.5, 1), name = 'Morphine vs. Control P < 0.05') + 
  my_theme + theme(legend.position = 'bottom') + 
  xlab('Human gene log2(fold-change)') +
  ylab('Monkey gene log2(fold-change)')
dev.off()



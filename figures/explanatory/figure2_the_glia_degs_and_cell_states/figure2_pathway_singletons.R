## packages for data table processing 
library(tidyverse)
library(stringr)
library(RColorBrewer)
library(data.table)
library(simplifyEnrichment)
library(ComplexHeatmap)
library(here)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

FIGDIR='figures/explanatory/figure2_the_glia_degs_and_cell_states'
dir.create(here(FIGDIR, 'plots', 'gsea'), recursive = T, showWarnings = F)
dir.create(here(FIGDIR, 'tables'), recursive = T, showWarnings = F)

###################################
# 0) pre-set colors and cell types 
celltypes = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome',
              'D1/D2-Hybrid', 'Interneuron',  'Astrocytes', 'Endothelial', 
              'Microglia', 'Mural', 'Oligos', 'Oligos_Pre')

neuron_types = c('D1-Matrix', 'D2-Matrix',  'D1-Striosome', 'D2-Striosome',
                 'D1/D2-Hybrid', 'Interneuron', 'Neuron') %>% make.names()

glia_types = c('Astrocytes', 'Endothelial', 'Microglia', 'Mural', 'Oligos', 
               'Oligos_Pre') %>% make.names()

in2mm<-25.4

###########################################
# 1) read in the pathway enrichment tables

## read in the clustered pathways, limit to unclustered pathways
singleton_df = here(FIGDIR, 'tables','figure2_clustered_glia_gsea_pathway_network.xlsx') %>% 
  readxl::read_xlsx() %>% filter(is.na(cluster_number)) %>% 
  filter(celltype %in% glia_types)

id_to_term = singleton_df %>% 
  mutate(description = gsub('[*]', '', description),
         description = gsub('cells', 'cell', description)) %>% 
  filter(!duplicated(pathway)) %>% dplyr::select(pathway, description) %>% deframe()

gsea_wide = singleton_df %>% arrange(celltype) %>% 
  pivot_wider(id_cols = pathway, names_from = 'celltype', values_from = 'NES', 
              values_fill = 0) %>% 
  column_to_rownames('pathway') %>% as.matrix()

go_id = rownames(gsea_wide)
go_term = id_to_term[go_id]

## create a singleton GO similarity matrix by gene overlaps
gene_list = gsea_df %>% group_by(pathway) %>% 
  summarise(value = paste(leadingEdge, collapse = ',')) %>% deframe() %>% 
  lapply(str_split,',') %>% lapply(unlist) %>% lapply(unique)
gene_list = gene_list[go_id]
mat = term_similarity(gene_list)

# hierarchicaly clustering of the pathway similarities
hr <- hclust(dist(mat, method = "euclidean"), method="complete") 
mycl.r <- cutree(hr, k=6)
mycl.r = setNames(letters[mycl.r], names(mycl.r))

## complex heatmap plot parameters
gp_title = gpar(fontsize = 6, fontface = 'bold'); gp_label = gpar(fontsize = 6)
u0 = unit(0, "mm"); u2 = unit(2, "mm"); u3 = unit(3, "mm")
pal <- RColorBrewer::brewer.pal(11, 'PiYG')

row_params = list(title_gp = gp_title, labels_gp = gp_label, 
                  column_gap = u2, grid_height = u2, grid_width = u3, 
                  legend_height = u3)

# create the word cloud row names of shared pathways
align_to = split(match(names(mycl.r), go_id), mycl.r)
term = split(go_term, mycl.r)
word_cloud = anno_word_cloud(align_to, term, max_words = 5, 
                             fontsize_range = c(2, 8),
                             bg_gp = gpar(fill = "#DDDDDD", col = "#CCCCCC"),
                             exclude_words = c('geneid', 'cancer'), 
                             link_width = unit(1, "mm"))


## make the plot
pdf(height = 70/in2mm, width = 59/in2mm)
Heatmap(gsea_wide, row_split= mycl.r, col = pal, show_heatmap_legend = F,
        show_row_names = F, show_column_names = T,
        cluster_rows = T, show_row_dend = F, 
        cluster_columns = T, show_column_dend = F, row_title = NULL,
        column_names_side = 'top', column_names_gp = gp_label,
        right_annotation = rowAnnotation(foo = word_cloud, annotation_width = unit(38, 'mm')))
dev.off()





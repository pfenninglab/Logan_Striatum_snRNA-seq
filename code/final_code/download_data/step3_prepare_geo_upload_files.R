library(tidyverse)
library(here)
library(Seurat)
library(SeuratDisk)

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }

#############################################
## 1) format the files for GEO submission
fq_fn = list.files(here('data/raw_data/fastq'), full.names = F, pattern = '.gz')

file_df = data.frame(path = fq_fn) %>% 
  mutate(Sample = basename(path) %>% ss('_merged'), 
         type = basename(path) %>% ss('_', 5)) %>% 
  pivot_wider(id_cols = Sample, names_from = type, values_from = path)

pd = readRDS(here('data/tidy_data/tables/LR_RM_OUD_snRNAseq_SampleInfo.rds')) %>% 
  dplyr::select(-c(X, Index.27:i5.index.seq, DSM.IV.OUD:DSM.IV.CUD)) %>% 
  rownames_to_column('Sample') %>% inner_join(file_df)

pd %>% writexl::write_xlsx(here('data/tidy_data/tables/LR_RM_OUD_snRNAseq_SampleInfo_fastq_files.xlsx'))

######################################################
## 2) annotate Seurat object for cellxgene submission
obj_merged = here('data/tidy_data/Seurat_projects', 
                  "BU_OUD_Striatum_refined_all_SeuratObj_N22.h5Seurat") %>% 
  LoadH5Seurat()

## set main counts for RNA
DefaultAssay(obj_merged) = 'RNA'
obj_merged@meta.data$celltype3 = 
  ifelse(grepl('Int', obj_merged$celltype3),
         "Interneuron", obj_merged$celltype3)
table(obj_merged$celltype3)

### add miscellaneous data
Misc(object = obj_merged, slot = "schema_version") = '3.0.0'
Misc(object = obj_merged, slot = "title") = 'Transcriptional responses of the human dorsal striatum in opioid use disorder implicates cell type-specifc programs'
Misc(object = obj_merged, slot = "batch_condition") = c('sex', 'tissue', 'disease', 'donor_id')

# pre-set cell types mapping
celltypes = c('D1-Matrix' = 'CL:4023026', # direct pathway medium spiny neuron
              'D2-Matrix' = 'CL:4023029', # indirect pathway medium spiny neuron
              'D1-Striosome' = 'CL:4023026',# direct pathway medium spiny neuron
              'D2-Striosome' = 'CL:4023029',# indirect pathway medium spiny neuron
              'D1/D2-Hybrid' = 'CL:1001474', # medium spiny neuron
              'Interneuron' = 'CL:0000498', # inhibitory interneuron
              'Astrocytes' = 'CL:0000127', # astrocyte
              'Endothelial' = 'CL:0000115',  # endothelial cell
              'Microglia' = 'CL:0000129', # microglial cell
              'Mural' = 'CL:0008034', # mural cell
              'Oligos' = 'CL:0000128', # oligodendrocyte
              'Oligos_Pre' = 'CL:0002453') #oligodendrocyte precursor cell

## create the mappings to required cell metadata
df = obj_merged[[]] %>% 
  mutate(
    # human NCBITaxon:9606
    organism_ontology_term_id = 'NCBITaxon:9606', organism = 'human', 
    # caudate UBERON:0001873, putamen UBERON:0001874
    tissue_ontology_term_id = ifelse(Region == 'Caudate','UBERON:0001873' ,'UBERON:0001874'), 
    tissue = Region, 
    # single nucleus RNA sequencing EFO:0009809
    assay_ontology_term_id = "EFO:0009922", assay = 'single nucleus RNA sequencing',
    # opioid dependence MONDO:0005530, normal PATO:0000461
    disease_ontology_term_id = ifelse(DSM.IV.OUD == 'OUD','MONDO:0005530' ,'PATO:0000461'), 
    disease = DSM.IV.OUD,
    # map the celltype ontology 
    cell_type_ontology_term_id = celltypes[celltype3], cell_type = celltype3,
    # African American HANCESTRO:0568, European/white HANCESTRO_0005
    self_reported_ethnicity_ontology_term_id = ifelse(Race == 'B', 'HANCESTRO:0568', 'HANCESTRO_0005'),
    self_reported_ethnicity = Race,
    # human adult stage
    development_stage_ontology_term_id = 'HsapDv:0000087', 
    development_stage = 'Adult', 
    # PATO:0000383 for female, PATO:0000384 for male, 
    sex_ontology_term_id = ifelse(Sex == 'F', 'PATO:0000383', 'PATO:0000384'),
    sex = Sex,
    # free text string, encouraged to be something not likely to be used in other studies
    donor_id = paste('LR_RM', Case, Sex, Race, Age, sep = '.'),
    suspension_type = 'nucleus', is_primary_data = T,
    )

## check the mappings look right
with(df, table(tissue_ontology_term_id, Region))
with(df, table(disease_ontology_term_id, DSM.IV.OUD))
with(df, table(cell_type_ontology_term_id, celltype3))
with(df, table(self_reported_ethnicity_ontology_term_id, Race))
with(df, table(sex_ontology_term_id, Sex))
with(df, table(donor_id, Case))

df2 = df %>% dplyr::select(-all_of(names(obj_merged[[]])))

## add the cell metadata to the object
obj_merged= AddMetaData(obj_merged, df2)
names(obj_merged[[]])
drop= c('Sex', 'Region',  'DSM.IV.OUD', 'celltype3', 'Race')
obj_merged@meta.data = obj_merged[[]] %>% 
  dplyr::select(-all_of(drop), -c(Index.27:integrated_snn_res.0.1))

## look at the gene metadata
genes_df = read_tsv('/home/bnphan/resources/genomes/GRCh38.p13/geneInfo.tab', 
                    skip = 1, col_names = c('ENSEMBL', 'feature_name', 'type')) 
genes_df = genes_df[match(rownames(obj_merged), genes_df$feature_name),] %>% 
  data.frame() %>% mutate(feature_is_filtered = F, 
                          feature_biotype = 'feature_name', 
                          feature_reference = "NCBITaxon:9606")
rownames(genes_df) = rownames(obj_merged)

## add the gene metadata to the features
obj_merged[["RNA"]] = 
  AddMetaData(obj_merged[["RNA"]], genes_df[rownames(obj_merged[["RNA"]]),])
obj_merged[["integrated"]] = 
  AddMetaData(obj_merged[["integrated"]], genes_df[rownames(obj_merged[["integrated"]]),])
obj_merged[["SCT"]] = 
  AddMetaData(obj_merged[["SCT"]], genes_df[rownames(obj_merged[["SCT"]]),])
newnames = genes_df %>% dplyr::select(feature_name, ENSEMBL) %>% deframe()

# RenameGenesSeurat  ------------------------------------------------------------------------------------
RenameGenesSeurat <- function(obj = ls.Seurat[[i]], newnames) { 
  # Replace gene names in different slots of a Seurat object. Run this before integration. Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.
  RNA <- obj@assays$RNA
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames[RNA@counts@Dimnames[[1]]]
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames[RNA@data@Dimnames[[1]]]
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames[RNA@scale.data@Dimnames[[1]]]
  obj@assays$RNA <- RNA

  SCT <- obj@assays$SCT
    if (length(SCT@counts)) SCT@counts@Dimnames[[1]] <- newnames[SCT@counts@Dimnames[[1]]]
    if (length(SCT@data)) SCT@data@Dimnames[[1]] <- newnames[SCT@data@Dimnames[[1]]]
    if (length(SCT@scale.data)) rownames(SCT@scale.data) <- newnames[rownames(SCT@scale.data)]
    if (length(SCT@var.features)) SCT@var.features <- newnames[SCT@var.features]
  obj@assays$SCT <- SCT

  integrated <- obj@assays$integrated
  if (length(integrated@var.features)) integrated@var.features <- newnames[integrated@var.features]
  if (length(integrated@counts)) rownames(integrated@counts) <- newnames[rownames(integrated@counts)]
  if (length(integrated@data)) rownames(integrated@data) <- newnames[rownames(integrated@data)]
  if (length(integrated@scale.data)) rownames(integrated@scale.data) <- newnames[rownames(integrated@scale.data)]
  obj@assays$integrated <- integrated
  return(obj)
}

obj_renamed = RenameGenesSeurat(obj = obj_merged, newnames) 
obj_renamed[['RNA']]
obj_renamed[['SCT']]
obj_renamed[['integrated']]

## save this verion to h5Seurat that can be converted to h5ad
out_fn =  here('data/tidy_data/geo_objects', 
               "BU_OUD_Striatum_refined_all_SeuratObj_N22.cellxgene.h5Seurat")
SaveH5Seurat(obj_renamed, filename = out_fn, overwrite = T)
Convert(out_fn, dest = "h5ad", overwrite = T)

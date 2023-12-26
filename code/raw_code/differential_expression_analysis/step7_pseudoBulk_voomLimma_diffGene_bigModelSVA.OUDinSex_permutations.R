## conda activate r4
## packages for data table processing 
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(optparse))

ss <- function(x, pattern, slot = 1, ...) { 
  sapply(strsplit(x = x, split = pattern, ...), '[', slot) }
options(stringsAsFactors = F)
options(repr.plot.width=11, repr.plot.height=8.5)

DATADIR='data/tidy_data/differential_expression_analysis'

###################################################
# 0) set the random seed from command line argument
option_list = list(
  make_option(c("-i", "--iter_num"), type="numeric", default= Sys.time() %>% as.numeric(), 
              help="number to set random seed", metavar="numeric")
)
opt = OptionParser(option_list=option_list) %>% parse_args();
print(paste('The random seed is:', opt$iter_num))
set.seed(opt$iter_num)

## check if this iter needs to be run, quit if already run
rdasDir =file.path(DATADIR, 'rdas', 'permutations2'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, paste0('OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.',opt$iter_num,'.rds'))
if(file.exists(save_res_fn)){
  print(paste('file exists, skipping:', save_res_fn))
  quit(save = 'no')
}

## main differential gene expression package
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(swfdr))

################################
# 1) load pseudobulk sce object
save_pseudobulk =here(DATADIR, 'rdas', 'BU_OUD_Striatum_refined_all_PseudoBulk_N22.sce2.rds')
pb = readRDS(save_pseudobulk)

####################################################
## 2) filter pseudobulk samples that have too few cells
pb = pb[, as.vector(pb$numCells > 15)]
pb = pb[, pb$celltype3 != 'Mural'] # drop mural cells b/c too few
pb$celltype3 = make.names(pb$celltype3) %>% as.factor()
pb$Region = as.factor(pb$Region)
pb$Sex = as.factor(pb$Sex)
pb$numCells = as.numeric(pb$numCells)

###############################################################################
## 3) permute the sex and OUD dx of the subjects and propagate across all cells
sex_mapping = colData(pb) %>% as.data.frame() %>% distinct(Case, Sex) %>% deframe()
sex_mapping = setNames(sample(unname(sex_mapping)), names(sex_mapping))
pb$Sex = sex_mapping[as.character(pb$Case)]

oud_mapping = colData(pb) %>% as.data.frame() %>% distinct(Case, DSM.IV.OUD) %>% deframe()
oud_mapping = setNames(sample(unname(oud_mapping)), names(oud_mapping))
pb$DSM.IV.OUD = oud_mapping[as.character(pb$Case)]

#################################################################
## 4) permute the sex of the subjects and propagate across all cells
## make an interaction term for all the combinations
pb$celltype_dx_rg_sex = interaction(pb$DSM.IV.OUD, pb$celltype3,pb$Sex, pb$Region) %>% 
  as.factor() %>% droplevels()

## interaction term w/o the DX
pb$celltype_rg_sex = interaction(pb$celltype3,pb$Sex, pb$Region) %>% 
  as.factor() %>% droplevels()

## construct design & contrast matrix regressing out the DetRate
design <- model.matrix(~ 0 + celltype_dx_rg_sex  + # term capturing the main effects
                         Age + PMI + RIN + numCells + cdr, # co-variates 
                       data = colData(pb))

## construct the null model, used in regressing out the various factors
design0 <- model.matrix(~ 0 +  celltype_rg_sex  + # term capturing the effects w/o Dx
                          Age + PMI + RIN  + numCells + cdr, # co-variates 
                        data = colData(pb))

####################################
# 5) normalization using voom-limma
y <- DGEList(counts = assays(pb)[[1]])

## filter out genes w/ low counts
A <- rowMeans(y$counts)
isexpr <- A > 5
y <- y[isexpr, , keep.lib.size = FALSE]

## filter out ribosomal genes, filter out mitochondria genes
keep.genes <- grep("^RP[SL]|^MT-",rownames(y), value = T, invert = T)
y = y[keep.genes, , keep.lib.size = FALSE]

# normalize counts
y <- calcNormFactors(y)

## voom precision weights and sample-wise quality weights normalization
v <- voomWithQualityWeights(y, design)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Case)

## recalculate weights after adjusting for correlated samples from same subject
v <- voomWithQualityWeights(y, design, block = colData(pb)$Case, 
                            correlation = cor$consensus)
cor <- duplicateCorrelation(v, design, block = colData(pb)$Case)

############################################################
# 6) Use surrogate variables to estimate unmodeled variation

## estimate the number of SVs from the adjusted
# (n.sv = num.sv(v$E, design, method="be", seed = set.seed(opt$iter_num))) #21
svobj = sva(v$E, design, design0, n.sv=21, B = 20)

## add the SVs to the model matrix
designSV = cbind(design, svobj$sv)
design0SV = cbind(design0, svobj$sv)

## recalculate sample quality weights after calculating the SVs
v <- voomWithQualityWeights(y, designSV, block = colData(pb)$Case, 
                            correlation = cor$consensus)
cor <- duplicateCorrelation(v, designSV, block = colData(pb)$Case)

## fit the model to get DEGs for various differences in the data
fit <- lmFit(v, designSV, block = colData(pb)$Case, correlation = cor$consensus)
fit <- eBayes(fit, robust = TRUE)

###########################################################
## 7) compute the differences b/t Dx within each cell type
celltypes=levels(factor(pb$celltype3 %>% make.names()))
designSV2 =designSV
colnames(designSV2) = make.names(colnames(designSV2))

## make the cell type contrasts in Males
con_celltypes_male = sapply(setNames(celltypes, celltypes),function(cell) {
  cell = colnames(designSV2) %>% make.names() %>% 
    str_subset(paste0('\\.', cell, '\\.')) %>% 
    str_subset(paste0('\\.M\\.')) ## this part here looks in just Males
  OUD = cell %>% str_subset('OUD') 
  CTL = cell %>% str_subset('CTL')
  
  N_SexF = OUD %>% length()
  OUD = OUD %>% paste(collapse = ' + ')
  
  N_SexM= CTL %>% length()
  CTL = CTL %>% paste(collapse = ' + ')
  
  paste('(',OUD,')/',N_SexF, '-(',CTL,')/',N_SexM)
})

## proportion of each cell type
df_prop = here('data/tidy_data/tables/BU_OUD_Striatum_refined_celltype3_proportions.txt') %>% 
  read_tsv(progress =F, show_col_types = F) %>% 
  deframe()
names(df_prop) = names(df_prop) %>% make.names()
df_prop = df_prop[celltypes]

ind_neur = grepl('^D|^Int',celltypes)
ind_glia = !grepl('^D|^Int',celltypes)

## create the contrasts for OUD effect Between all cells or major classes
con_groups_male = c('All' = paste0('(', con_celltypes_male,')*', df_prop) %>% paste(collapse = ' + '), 
               'Neuron' = paste0('(', con_celltypes_male[ind_neur],')*', 
                                 df_prop[ind_neur]/sum(df_prop[ind_neur])) %>% paste(collapse = ' + '), 
               'Glia' =  paste0('(', con_celltypes_male[ind_glia],')*', 
                                df_prop[ind_glia]/sum(df_prop[ind_glia])) %>% paste(collapse = ' + '))

## swap the male-OUD effects w/ contracts for female-OUD effects
con_male = c(con_groups_male, con_celltypes_male)
con_female = con_male %>% str_replace_all('\\.M\\.', '.F.') %>% 
  setNames(names(con_male) %>% paste0('#SexF'))
names(con_male) = paste0(names(con_male), '#SexM')
cont.matrix <- makeContrasts(contrasts= c(con_female, con_male), levels=designSV2)
rownames(cont.matrix) = colnames(designSV)

## refit the model based on these contrasts
fit2 <- contrasts.fit(fit, cont.matrix) %>% eBayes()

## compute the DEGs from these contrasts
deg_list = lapply(setNames(colnames(cont.matrix), names(c(con_female, con_male))), 
                  function(coef){
                    topTable(coef = coef, fit =fit2, n=Inf) %>% arrange(P.Value) %>% 
                      ## use SWFDR to increase power of detecting DEGs based on avg expression covariate
                      ## https://pubmed.ncbi.nlm.nih.gov/30581661/
                      mutate(adj.P.Val.Within =  lm_qvalue(P.Value, X=AveExpr)$q) %>%
                      dplyr::select(-adj.P.Val) %>% rownames_to_column('gene') 
                  })

## FDR correction Between all tests
deg_list = deg_list %>% data.table::rbindlist(idcol = 'group') %>% 
  mutate(adj.P.Val.Between =  lm_qvalue(P.Value, X=AveExpr)$q, 
         celltype = ss(group, '#', 1), Sex = ss(group, '#', 2)) %>%
  dplyr::select(-group) %>% pivot_longer(-c(gene, celltype, Sex)) %>% 
  pivot_wider(id_cols = c(celltype, gene), names_from = c(name, Sex), 
              values_from = value) %>% 
  split(f = .$celltype)

## calculate the interaction score between the 2 conditions
deg_list = deg_list %>% lapply(function(x){
  x %>% 
    ## calculate a number that gives both statistical signif + effect size
    mutate(dir_SexF = -log10(P.Value_SexF) * sign(logFC_SexF),
           dir_SexM = -log10(P.Value_SexM) * sign(logFC_SexM),
           ## this will score genes most different b/t 2 conditions
           dir_difference = (dir_SexF - dir_SexM)) %>% 
    arrange(desc(abs(dir_difference)))
}) %>% bind_rows() %>% 
  ## skinny version of the dataframe
  dplyr::select(celltype, gene, starts_with('dir'), everything()) %>%
  ## add the permutation ID to the results
  mutate(iter_num = opt$iter_num)

####################################################################
## 8) save the output of voom_limma differential state analyses
saveRDS(deg_list, save_res_fn)
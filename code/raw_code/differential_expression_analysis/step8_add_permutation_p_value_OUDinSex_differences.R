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
# 1) read in the main results from the sex-specific results
rdasDir =file.path(DATADIR, 'rdas'); dir.create(rdasDir, showWarnings = F)
save_res_fn = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.rds')

deg_df = readRDS(save_res_fn) %>% bind_rows() %>% 
  mutate(dir_num_extreme = 0, # number of observations with a more extreme value than the observed
         dir_num_perm = 0, # number of permutations
         tmp = paste(celltype,gene, sep='_'))

## get the unpermuted absolute values
abs_val1 = deg_df %>% 
  mutate(tmp = paste(celltype,gene, sep='_')) %>% 
  dplyr::select(tmp, dir_difference) %>% deframe() %>% abs()

#################################################################
# 2) read in each permutation and update the deg_df with extreme values
files =here(DATADIR, 'rdas', 'permutations2') %>% list.files(full.names = T)

# read in ncores worth of files at a time
ncore = 8
files_list = split(files, ceiling(seq_along(files)/ncore))

for(i in seq_along(files_list)){
  startTime <- Sys.time() 
  
  print(paste(i, 'of', length(files_list)))
  
  ## read in parallel chunks of permuted data
  extreme_vec= parallel::mclapply(files_list[[i]], function(file){
    abs_val2 = readRDS(file) %>% 
      mutate(tmp = paste(celltype,gene, sep='_')) %>% 
      dplyr::select(tmp, dir_difference) %>% deframe() %>% abs()
    
    ## re-order to match the order of the observed data
    abs_val2 = abs_val2[names(abs_val1)]
    
    ## count the number of permuted values that are more extreme than the observed
    return(abs_val2 > abs_val1)
  }, mc.cores = ncore) %>% bind_cols() %>% 
    apply(1, sum)
  
  deg_df = deg_df %>% 
    mutate(
      ## keep running track how many times permuted values are more extreme than the observed
      dir_num_extreme = dir_num_extreme + extreme_vec, 
      ## keep running track of how many permutations have been done
      dir_num_perm = dir_num_perm + length(files_list[[i]]), 
    )
  
  # prints recorded time per iteration
  endTime <- Sys.time() 
  print(endTime - startTime)
}

#################################################################
# 3) calculate the permutation p-value for each gene and cell type
diff_list = deg_df %>% 
  mutate(dir_difference_permP = dir_num_extreme/dir_num_perm) %>% 
  mutate(dir_difference_permP = pmax(1/dir_num_perm, dir_difference_permP)) %>% 
  dplyr::select(-tmp) %>% 
  split(f = .$celltype)

save_res_fn2 = here(rdasDir, 'OUD_Striatum_voom_limma_bigModelSVA_N22.OUDwinSex.permP.rds')
saveRDS(diff_list, save_res_fn2)


## see how many genes are significant at different p-value thresholds
sapply(diff_list, function(df) df %>% filter(dir_difference_permP < 0.05) %>% nrow())
sapply(diff_list, function(df) df %>% filter(dir_difference_permP < 0.01) %>% nrow())
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import spmatrix,coo_matrix,csc_matrix
from Bio import Phylo
from io import StringIO
import sys, os
sys.path.append("/home/bnphan/src/TMNT/tmnt_algorithm")
import TMNT

DATADIR="data/tidy_data/hierarchical_deg_analyses/"
save_file = "data/tidy_data/Seurat_projects/OUD_Striatum_refined_all_SeuratObj_N22.h5ad"
striatum_ann = sc.read_h5ad(save_file)

## check has the cell type & Dx
striatum_ann.obs['level1'] = striatum_ann.obs['level1'].astype('string')
striatum_ann.obs.groupby(['level1']).size()

striatum_ann.obs['level2'] = striatum_ann.obs['level2'].astype('string')
striatum_ann.obs.groupby(['level2']).size()

## the cell type, Dx tree, and 
cell_type_tree = Phylo.read(DATADIR + '/tables/'+ 'OUD_Striatum_refined_celltype.byComp.nh', "newick")
oud_dx_tree = Phylo.read(StringIO("(OUD)CTL;"), "newick") ## 2 cases
input_trees = [oud_dx_tree,cell_type_tree]

## get the top variable genes
select_gene_df = striatum_ann.var.loc[striatum_ann.var["highly_variable"].values,:]
select_gene_df.to_csv(DATADIR + "/tables/OUD_Striatum_refined_scanpy_N22_selected_genes.csv")
selected_striatum_matrix = coo_matrix(striatum_ann.X[:,striatum_ann.var["highly_variable"].values])

## run TMNT w/ different combinations of lasso models & intercepts
lasso_penalty_list = ["LarsAIC","LarsBIC","Lasso.01","ElasticNet","Lasso.001","LarsCV","LassoCV","Lasso.1"]
lasso_penalty_list = ["LarsAIC"]
intercept_list= [True,False]
intercept_list= [True]

print('getting started')
lasso_penalty = "LarsAIC"
intercept = True

run_id = lasso_penalty+"_PCA-"+str(intercept)
fname = DATADIR + "/tmnt_runs/"+run_id + '/'
os.mkdir(fname)
print(run_id)
mc_TMNT = TMNT.NestedTree(selected_striatum_matrix, select_gene_df, 
    striatum_ann.obs, input_trees, model_name=lasso_penalty, intercept=intercept)
mc_TMNT.ancestor_matrix.to_csv(DATADIR + "/tables/tmnt_ancestor_df.csv")
mc_TMNT.em_algorithm_checkpoints(epsilon = 1e-10, save_folder = fname)
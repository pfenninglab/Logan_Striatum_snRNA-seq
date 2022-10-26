import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
import sys

DATADIR="data/tidy_data/hierarchical_deg_analyses/"
save_file = "data/tidy_data/Seurat_projects/OUD_Striatum_refined_all_SeuratObj_N22.h5ad"
striatum_ann = sc.read_h5ad(save_file)

## relabel the columns for multi-level hierarchy
dx_dict = {0 : "CTL", 1:"OUD"} # b/c Seurat to h5ad conversion doesn't like factors
striatum_ann.obs["level1"] = pd.Series([dx_dict[region_name] for region_name in striatum_ann.obs["DSM.IV.OUD"].values ],dtype=str,index=striatum_ann.obs.index)
striatum_ann.obs["level2"] = pd.Series(striatum_ann.obs["celltype3"])

## check lables is correct
striatum_ann.obs.groupby(["level1", "Dur.OUD"]).size()
striatum_ann.obs.groupby("level2").size()

## scanpy preprocessing
striatum_ann.var_names_make_unique()  
sc.pp.filter_cells(striatum_ann, min_genes=200)
sc.pp.filter_genes(striatum_ann, min_cells=3)
sc.pp.normalize_total(striatum_ann, target_sum=1e4)
sc.pp.log1p(striatum_ann)
sc.pp.highly_variable_genes(striatum_ann,n_top_genes=16000)

## Save the pre-processed h5ad file
striatum_ann.__dict__['_raw'].__dict__['_var'] = striatum_ann.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
striatum_ann.write(save_file)
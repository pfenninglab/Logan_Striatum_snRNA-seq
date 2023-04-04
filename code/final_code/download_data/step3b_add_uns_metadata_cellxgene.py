import numpy as np
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
print(ad.__version__)

h5_fn = 'data/tidy_data/geo_objects/BU_OUD_Striatum_refined_all_SeuratObj_N22.cellxgene.h5ad'
adata = ad.read(h5_fn)
adata.var_names_make_unique()

adata.uns["schema_version"] = '3.0.0'
adata.uns["title"] = 'Transcriptional responses of the human dorsal striatum in opioid use disorder implicates cell type-specifc programs'
adata.uns["batch_condition"] = ['sex', 'tissue', 'disease', 'donor_id']

adata.write(h5_fn, compression = 'gzip')

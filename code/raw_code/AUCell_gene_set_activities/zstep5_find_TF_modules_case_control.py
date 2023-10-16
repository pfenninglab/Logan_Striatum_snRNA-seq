import os
import glob
import pickle
import pandas as pd
import numpy as np
import loompy as lp

from dask.diagnostics import ProgressBar

from arboreto.utils import load_tf_names

from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies, load_motifs
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell
import pyscenic

import seaborn as sns

def name(fname):
    return os.path.splitext(os.path.basename(fname))[0]

DATA_FOLDER="/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/tidy_data/AUCell_gene_set_activities"
DATABASE_FOLDER = "/home/bnphan/resources/SCENIC"
SCHEDULER="123.122.8.24:8786"
DATABASES_GLOB = os.path.join(DATABASE_FOLDER, "hg38_*full_tx_v10_clust.genes_vs_motifs.rankings.feather")
MOTIF_ANNOTATIONS_FNAME = "/home/bnphan/resources/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"
MM_TFS_FNAME = "/home/bnphan/resources/SCENIC/allTFs_hg38.txt"
SC_EXP_FNAME = os.path.join(DATA_FOLDER, 'regulon', "OUD_Striatum_refined_N22.loom")
REGULONS_FNAME = os.path.join(DATA_FOLDER, 'regulon', "OUD_Striatum_refined_N22.regulons.pickle")
ADJACENCIES =  os.path.join(DATA_FOLDER, 'grnboost2', "Combined_OUD_Striatum.GRNBoost2.tsv")

#################################################################
## 1) grab the expression matr
ds = lp.connect(SC_EXP_FNAME)
ex_matrix = pd.DataFrame(data=ds[:, :], index=ds.ca.CellID, columns=ds.ra.Gene)  # 1st row as the column names
ex_matrix =ex_matrix.T
ds.close()

tf_names = load_tf_names(MM_TFS_FNAME)

## import the TF-target adjacency table and create modules
adjacencies = pd.read_table(ADJACENCIES)
modules = list(modules_from_adjacencies(adjacencies, ex_matrix))



db_fnames = glob.glob(DATABASES_GLOB)
dbs = [RankingDatabase(fname=fname, name=name(fname)) for fname in db_fnames]
a.geneset

file = open('/home/bnphan/resources/SCENIC/hg38_refseq80_genes.txt','w')
for item in a.geneset:
    file.write(item+"\n")

file.close()

#################################################################
## Calculate a list of enriched motifs and the corresponding target genes for all modules.
with ProgressBar():
    df = prune2df(dbs, modules, MOTIF_ANNOTATIONS_FNAME)


######################################################
## Create regulons from this table of enriched motifs.
regulons = df2regulons(df)

# Save the enriched motifs and the discovered regulons to disk.
df.to_csv(MOTIFS_FNAME)
with open(REGULONS_FNAME, "wb") as f:
    pickle.dump(regulons, f)



auc_mtx = aucell(ex_matrix, regulons, num_workers=4)

AUCELL=os.path.join(DATA_FOLDER, 'regulon', "OUD_Striatum_refined_N22.aucell.tsv")

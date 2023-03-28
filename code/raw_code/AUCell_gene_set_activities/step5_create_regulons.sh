#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen_bigmem
#SBATCH --time 6:00:00
##SBATCH -w compute-1-40
#SBATCH --job-name=ctx
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=125G
#SBATCH --dependency=afterok:3315375
#SBATCH --error=logs/pyscenic_cistarget_%A.txt
#SBATCH --output=logs/pyscenic_cistarget_%A.txt

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/AUCell_gene_set_activities
TF_LIST=/home/bnphan/resources/SCENIC/allTFs_hg38.txt
DBNAME1=/home/bnphan/resources/SCENIC/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
DBNAME2=/home/bnphan/resources/SCENIC/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather
ANNOT=/home/bnphan/resources/SCENIC/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
mkdir -p $DATADIR/grnboost2 $DATADIR/regulon
cd $PROJDIR/code/raw_code/AUCell_gene_set_activities

source ~/.bashrc; mamba activate pyscenic

###################################################
## 1) prune the TF-gene targets and create regulons
MODULE=$DATADIR/grnboost2/Combined_OUD_Striatum.GRNBoost2_filtered.tsv
EXPRMAT=$DATADIR/regulon/OUD_Striatum_refined_N22.loom
REGULON=$DATADIR/regulon/OUD_Striatum_refined_N22.regulons.tsv

pyscenic ctx -o ${REGULON} \
--mode custom_multiprocessing \
--expression_mtx_fname ${EXPRMAT} \
--annotations_fname ${ANNOT} \
--cell_id_attribute "CellID" \
--gene_attribute "Gene" \
--thresholds 0.90 \
--num_workers 8 \
${MODULE} ${DBNAME1} ${DBNAME2}
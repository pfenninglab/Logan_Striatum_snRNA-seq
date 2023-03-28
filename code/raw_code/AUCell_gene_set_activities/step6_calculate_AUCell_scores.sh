#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --time 3:00:00
#SBATCH --job-name=aucell
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=45G
#SBATCH --error=logs/pyscenic_aucell_%A.txt
#SBATCH --output=logs/pyscenic_aucell_%A.txt

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/AUCell_gene_set_activities
mkdir -p $DATADIR/grnboost2 $DATADIR/regulon
cd $PROJDIR/code/raw_code/AUCell_gene_set_activities

source ~/.bashrc; mamba activate pyscenic

###############################################################
## 1) calculate the single cell AUCell scores across regulons
EXPRMAT=$DATADIR/regulon/OUD_Striatum_refined_N22.loom
REGULON=$DATADIR/regulon/OUD_Striatum_refined_N22.regulons.tsv
AUCELL=$DATADIR/regulon/OUD_Striatum_refined_N22.aucell.tsv

pyscenic aucell -o ${AUCELL} \
--cell_id_attribute "CellID" \
--gene_attribute "Gene" \
--num_workers 16 \
${EXPRMAT} ${REGULON}

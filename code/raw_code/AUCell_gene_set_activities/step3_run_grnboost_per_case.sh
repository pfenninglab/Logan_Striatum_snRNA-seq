#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --exclude=compute-1-11,compute-1-12,compute-1-35
##SBATCH --time 3-00:00:00
#SBATCH --job-name=grnboost
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=46G
#SBATCH --error=logs/pyscenic_grnboost_%A_%a.txt
#SBATCH --output=logs/pyscenic_grnboost_%A_%a.txt
##SBATCH --array=8-13,20-25,32-37,44-49,56-61,68-73%4
#SBATCH --array=0-119%6

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/AUCell_gene_set_activities
TF_LIST=/home/bnphan/resources/SCENIC/allTFs_hg38.txt
mkdir -p $DATADIR/grnboost2
cd $PROJDIR/code/raw_code/AUCell_gene_set_activities

source ~/.bashrc; mamba activate pyscenic; 

########################################################
## 1) run the GRNBoost2 w/ pysenic N times for this file
N=`expr $SLURM_ARRAY_TASK_ID % 12 + 1`
ITER=`expr $SLURM_ARRAY_TASK_ID / 12 + 1`
EXPRMAT=$(ls $DATADIR/loom/* | sed "${N}q;d")
LABEL=$(basename $EXPRMAT .loom | sed 's/.*\.//g')

OUTFILE=$DATADIR/grnboost2/OUD_Striatum_GRNBoost2.${LABEL}.${ITER}.tsv
echo "Working on Sample ID:${LABEL} for the ${ITER}th run."
if [[ ! -f $OUTFILE || $(cat $OUTFILE) == '' ]]; then
arboreto_with_multiprocessing.py \
-o ${OUTFILE} -m grnboost2 \
--seed ${SLURM_ARRAY_TASK_ID} \
--cell_id_attribute "CellID" \
--gene_attribute "Gene" \
--num_workers 15 \
${EXPRMAT} ${TF_LIST}
fi
exit

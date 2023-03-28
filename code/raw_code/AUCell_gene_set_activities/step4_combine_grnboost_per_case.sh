#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 3-00:00:00
#SBATCH --job-name=combine
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/combine_grnboost_%A_%a.txt
#SBATCH --output=logs/combine_grnboost_%A_%a.txt
#SBATCH --array=1-12

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/AUCell_gene_set_activities
TF_LIST=/home/bnphan/resources/SCENIC/allTFs_hg38.txt
mkdir -p $DATADIR/grnboost2
cd $PROJDIR/code/raw_code/AUCell_gene_set_activities


########################################################
## 1) find replicate runs of GRNboost for a given sample
EXPRMAT=$(ls $DATADIR/loom/* | sed "${SLURM_ARRAY_TASK_ID}q;d")
LABEL=$(basename $EXPRMAT .loom | sed 's/.*\.//g')
REP=$(ls $DATADIR/grnboost2/OUD_Striatum_GRNBoost2.${LABEL}.*.tsv| wc -l|cut -f1)

## find the TF-target gene sample w/ more than N times ()
## average the importance scores in 3rd column
COUNT=$DATADIR/grnboost2/Avg_GRNBoost2.${LABEL}.tsv
echo -e "TF\ttarget\timportance" > $COUNT
awk FNR!=1 $DATADIR/grnboost2/OUD_Striatum_GRNBoost2.${LABEL}.*.tsv | \
awk -F '\t' -v REP=$M '
BEGIN { OFS = "\t"} 
NR>1{ $4 =$1"#"$2; arr[$4] += $3; count[$4] += 1} 
END{ for (a in arr) { if (count[a] >= REP * 0.8){
print a OFS arr[a]/count[a]
}}}' | tr '#' '\t' | sort -k1 -k2 >> $COUNT 
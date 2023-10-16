#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 3-00:00:00
#SBATCH --job-name=permute
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --error=logs/permut_sex_in_oud_degs_%A_%a.txt
#SBATCH --output=logs/permut_sex_in_oud_degs_%A_%a.txt
#SBATCH --array=1-1000

source ~/.bashrc
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
CODEDIR=$PROJDIR/code/raw_code/differential_expression_analysis/
# cd $CODEDIR

## update the array

for ROUND in  $(seq 0 1000 50000);
do
ID=`expr $SLURM_ARRAY_TASK_ID + $ROUND`

## run one iteration of the permutation
Rscript \
--vanilla $CODEDIR/step7_pseudoBulk_voomLimma_diffGene_bigModelSVA.OUDinSex_permutations.R \
--iter_num $ID
done

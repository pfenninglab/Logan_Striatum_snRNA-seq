#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen_bigmem
#SBATCH --time 3-00:00:00
#SBATCH --job-name=combine
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --dependency=afterok:3315363
#SBATCH --mem=48G
#SBATCH --error=logs/combine_all_%A.txt
#SBATCH --output=logs/combine_all_%A.txt

###################
## 1) set up paths 
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/AUCell_gene_set_activities
TF_LIST=/home/bnphan/resources/SCENIC/allTFs_hg38.txt
mkdir -p $DATADIR/grnboost2
cd $PROJDIR/code/raw_code/AUCell_gene_set_activities

##################################################
## merge the TF-target importances across samples
REP=12 ## number of samples
COUNT=$DATADIR/grnboost2/Combined_OUD_Striatum.GRNBoost2.tsv
echo -e "TF\ttarget\timportance" > $COUNT
awk FNR!=1 $DATADIR/grnboost2/Avg_GRNBoost2.*.tsv | \
awk -F '\t' '
BEGIN { OFS = "\t"} 
NR>1{ $4 =$1"#"$2; arr[$4] += $3; count[$4] += 1} 
END{for (a in arr) { if (count[a] >= REP * 0.8){
print a OFS arr[a]/count[a]
}}}' | tr '#' '\t' >> $COUNT


##########################################################################
## filter out TF-target pairings for genes in the SCENIC ranking databases
COUNT2=$DATADIR/grnboost2/Combined_OUD_Striatum.GRNBoost2_filtered.tsv
GENES=/home/bnphan/resources/SCENIC/hg38_refseq80_genes.txt
echo -e "TF\ttarget\timportance" > $COUNT2
awk 'NR==FNR{_[$1];next} ($2 in _)' $GENES $COUNT >> $COUNT2

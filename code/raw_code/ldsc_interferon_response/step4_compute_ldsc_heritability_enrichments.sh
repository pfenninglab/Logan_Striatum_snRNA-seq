#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time=3-0:00:00
#SBATCH --job-name=ldsc
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/est_herit_%A_%a.txt
#SBATCH --output=logs/est_herit_%A_%a.txt
#SBATCH --array=1-78

log2results() {
awk -F '\t' '/Total Observed scale h2*/{flag=1;next}/Lambda GC/{flag=0}flag' ${OUT}.log | \
# append line with leading white space to previous
sed ':a;N;$!ba;s/\n //g'| \
# parse Partitioned heritability rownames, starts w/ colon
awk -F":" -v OFS='\t' '{gsub("[[:space:]]", "_", $1);gsub(":", "", $0); print}' | \
# transpose row to columns
awk -v OFS='\t' '
{ 
for (i=1; i<=NF; i++){
a[NR,i] = $i
}
}
NF>p { p = NF }
END {    
for(j=1; j<=p; j++){
str=a[1,j]
for(i=2; i<=NR; i++){
str=str"\t"a[i,j];
}
print str
}
}' | awk -v OFS='\t' '{gsub("_[0-9]+$", "", $1); gsub("L2$", "", $1); print}' | \
gzip > ${OUT}.results.gz
}

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
SNPLIST=${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist
CODEDIR=$SETDIR/code/raw_code/ldsc_interferon_response
DATADIR=$SETDIR/data/tidy_data/ldsc_interferon_response
cd $CODEDIR; 
source activate ldsc

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' ${GWASDIR}/gwas_list_sumstats.tsv)
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $2}' ${GWASDIR}/gwas_list_sumstats.tsv)
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' ${GWASDIR}/gwas_list_sumstats.tsv)
OUTDIR=${DATADIR}/prop_herit; mkdir -p $OUTDIR

###################################################################
# run LD score regression over the Corces caudate binary annotations
CTS_FN=${DATADIR}/IFN_Macrophages_Monocytes.${POP}_hg38.ldcts
for IND in $(wc -l $CTS_FN | cut -d ' ' -f1 | xargs seq ); do
LAB=$(awk -v IND=$IND 'FNR == IND {print $1}' $CTS_FN)
if [[ $LAB != '' ]]; then
CELLTYPES=$(awk -v IND=$IND 'FNR == IND {print $2}' $CTS_FN)
CELL2=$(echo $CELLTYPES | cut -d ',' -f1)
OUT=${OUTDIR}/IFN_Macrophages_Monocytes.${GWAS_Label}.${POP}.${LAB}
if [[ ! -f "${OUT}.results.gz" || "${OUT}.results.gz" -ot ${CELL2}.1.l2.M ]]; then
ldsc.py --h2 $GWAS --w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--ref-ld-chr ${CELLTYPES},${GWASDIR}/1000G_ALL_Phase3_hg38_files/baselineLD_v2.2/baselineLD_v2.2.${POP}. \
--print-coefficients --out ${OUT} 
log2results; rm ${OUT}.log
if [[ $(zcat ${OUT}.results.gz) == '' ]]; then rm ${OUT}.results.gz; fi
fi; fi; done

## extract only the top cell type enrichment
OUTFILE=${OUTDIR}/IFN_Macrophages_Monocytes.${GWAS_Label}.${POP}.agg.gz
if [[ ! -f $OUTFILE || \
	$(zcat ${OUTFILE}) == '' || \
	"$OUTFILE" -ot ${CELL2}.1.l2.M ]]; then
gunzip ${OUTDIR}/IFN_Macrophages_Monocytes.${GWAS_Label}.${POP}*.results.gz
awk ' FNR == 2 || NR==1 {print}' \
${OUTDIR}/IFN_Macrophages_Monocytes.${GWAS_Label}.${POP}*.results |
gzip > ${OUTFILE}
gzip ${OUTDIR}/IFN_Macrophages_Monocytes.${GWAS_Label}.${POP}*.results
fi




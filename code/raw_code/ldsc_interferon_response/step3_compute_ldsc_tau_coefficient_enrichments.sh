#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1
#SBATCH --time 3-00:00:00
#SBATCH --job-name=ldsc_ifn_tau
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --error=logs/ldsc_ifn_tau_%A_%a.txt
#SBATCH --output=logs/ldsc_ifn_tau_%A_%a.txt
#SBATCH --array=1-78

# get the GWAS for this array job
SETDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
SNPLIST=${GWASDIR}/1000G_EUR_Phase3_GRCh38_files/hapmap3_snps/w_hm3.noMHC.snplist
CODEDIR=$SETDIR/code/raw_code/ldsc_interferon_response
DATADIR=$SETDIR/data/tidy_data/ldsc_interferon_response
cd $CODEDIR; source activate ldsc

# get the GWAS and reference population
GWAS=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $1}' /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas_list_sumstats.tsv )
POP=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND + 1) {print $2}' /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas_list_sumstats.tsv )
GWAS_Label=$(awk -F '\t' -v IND=${SLURM_ARRAY_TASK_ID} 'NR==(IND+ 1) {print $3}' /projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments/gwas_list_sumstats.tsv )

OUTDIR=${DATADIR}/enrichments_naive_bg; mkdir -p $OUTDIR

#################################################################################
# run LD score regression over the Stauffer striatum cell type binary annotations
CTS_FN1=${DATADIR}/IFN_Macrophages_Monocytes_Naive_BG.${POP}_hg38.ldcts
if [[ ! -f "$OUTDIR/IFN_Macrophages_Monocytes_Naive_BG.${GWAS_Label}.${POP}.cell_type_results.txt" ]]; then
ldsc.py --ref-ld-chr-cts $CTS_FN1 \
--ref-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/baseline_v1.1/baseline_v1.1.${POP}. \
--w-ld-chr ${GWASDIR}/1000G_ALL_Phase3_hg38_files/weights/1000G.${POP}.weights.hm3_noMHC. \
--h2-cts $GWAS --out $OUTDIR/IFN_Macrophages_Monocytes_Naive_BG.${GWAS_Label}.${POP}
fi

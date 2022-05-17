#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --time 3-00:00:00
#SBATCH --job-name=download
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=25G
#SBATCH --error=logs/download_mkfastq_%A_%a.txt
#SBATCH --output=logs/download_mkfastq_%A_%a.txt
#SBATCH --array=7%1

PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
cd $PROJDIR

## make sure these illumina download IDs are right
ILLUMINA_ID=(215462301 232949763 233627432 234243031 234543341 234860637 235352130 235650445)
RUNS=("Pilot" "Run1" "Run2" "Run3" "Run4" "Run5" "Run6" "Run7")
RUN=${RUNS[${SLURM_ARRAY_TASK_ID}]}

mkdir -p ${PROJDIR}/data/raw_data/bcl 

### download BCL files, grab run IDs from base space account w/ Ryan Logan's login
if [[ ! -f ${PROJDIR}/data/raw_data/bcl/BU_snRNAseq_${RUN}/RTAComplete.txt ]]; then
echo "Downloading the BCL files of the sequencing runs."
bs download run -i ${ILLUMINA_ID[${SLURM_ARRAY_TASK_ID}]} \
-o ${PROJDIR}/data/raw_data/bcl/BU_snRNAseq_${RUN}
fi


### make fastqs from the BCL sequencing runs
# if [[ ! -d ${PROJDIR}/data/raw_data/fastq/$RUN ]]; then
echo 'Making fastq files'
cellranger mkfastq --delete-undetermined --id=logs_Run_${RUN} \
--run=${PROJDIR}/data/raw_data/bcl/BU_snRNAseq_${RUN} \
--csv=${PROJDIR}/data/tidy_data/tables/LR_RM_BUMSR_Library_Index_Key_Run1_Run2.csv \
--project=${RUN} --output-dir ${PROJDIR}/data/raw_data/fastq
else
echo "Fastq folder already made at ${PROJDIR}/data/raw_data/fastq/${RUN}."
# fi
#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pool1,pfen1
#SBATCH --time 1-00:00:00
#SBATCH --job-name=clean_up
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=45G
#SBATCH --error=logs/STARsolo_%A_%a.txt
#SBATCH --output=logs/STARsolo_%A_%a.txt
#SBATCH --array=1-24%4
source ~/.bashrc

###########################################
# get sample name in demultiplex csv table
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=${PROJDIR}/data/raw_data
FASTQDIR=${DATADIR}/fastq
FASTQDIR2=${DATADIR}/fastq_merged
TMPDIR=/scratch/${USER}

mkdir -p $TMPDIR $FASTQDIR2

##########################
# get the sample file name
SAMPLE_ID=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $2}' ${PROJDIR}/data/tidy_data/tables/LR_RM_BUMSR_Library_Index_Key_Run1_Run2.csv )

# check if merged fastq files already made
if [[ -f $FASTQDIR2/${SAMPLE_ID}_merged_R2_001.fastq.gz ]]; then exit; fi 

cd $TMPDIR
rm $TMPDIR/${SAMPLE_ID}*

###########################################################
# copy over the fastq files, preserving Run file structure
for DIR in `ls ${FASTQDIR}/`; do
echo "Looking in ${DIR}"; mkdir -p ${DIR}
rsync -Paq ${FASTQDIR}/${DIR}/${SAMPLE_ID}_*_001.fastq.gz $TMPDIR/$DIR
done

###############################################################
# merge the fastq files across runs and lanes into one big file
I1_MERGED=${SAMPLE_ID}_merged_I1_001.fastq.gz 
ls */${SAMPLE_ID}_*I1_001.fastq.gz
cat */${SAMPLE_ID}_*I1_001.fastq.gz > $I1_MERGED

R1_MERGED=${SAMPLE_ID}_merged_R1_001.fastq.gz 
ls */${SAMPLE_ID}_*R1_001.fastq.gz
cat */${SAMPLE_ID}_*R1_001.fastq.gz > $R1_MERGED

R2_MERGED=${SAMPLE_ID}_merged_R2_001.fastq.gz 
ls */${SAMPLE_ID}_*R2_001.fastq.gz
cat */${SAMPLE_ID}_*R2_001.fastq.gz > $R2_MERGED

#####################################
## move back to main storage location
rsync -Paq --remove-source-files ${SAMPLE_ID}_*.fastq.gz $FASTQDIR2


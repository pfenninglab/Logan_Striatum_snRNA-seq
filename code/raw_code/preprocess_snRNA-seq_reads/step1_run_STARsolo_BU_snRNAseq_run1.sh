#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --time 3-00:00:00
#SBATCH --job-name=STARsolo
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=80G
#SBATCH --error=logs/STARsolo_%A_%a.txt
#SBATCH --output=logs/STARsolo_%A_%a.txt
#SBATCH --array=1-24%4

source ~/.bashrc

###########################################
# get sample name in demultiplex csv table
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=${PROJDIR}/data/raw_data
FASTQDIR=/data/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/raw_data
SOLODIR=$DATADIR/STARsolo_out
TMPDIR=/scratch/${USER}

mkdir -p $TMPDIR $SOLODIR; cd $SOLODIR

# perform mapping & read quantification
GENOME_DIR=/home/bnphan/resources/genomes/GRCh38.p13
BARCODES=/home/bnphan/src/cellranger-6.0.0/lib/python/cellranger/barcodes/3M-february-2018.txt

# sample-specific files
SAMPLE_ID=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $2}' \
${PROJDIR}/data/tidy_data/tables/LR_RM_BUMSR_Library_Index_Key_Run1_Run2.csv )

##################################################
# cDNA fragment in Read2, cell barcode in Read1

cd $TMPDIR
rm $TMPDIR/${SAMPLE_ID}*

if [[ ! -f "$PROJDIR/data/raw_data/STARsolo_out/${SAMPLE_ID}.Log.final.out" ]]; then

# copy over the fastq files, preserving Run file structure
for DIR in `ls ${FASTQDIR}/`; do
echo "Looking in ${DIR}"; mkdir -p ${DIR}
rsync -Paq ${FASTQDIR}/${DIR}/${SAMPLE_ID}_*R?_001.fastq.gz $TMPDIR/$DIR
done

# find all the files
cDNA_FASTQ=$(ls */${SAMPLE_ID}_*R2_001.fastq.gz | tr '[[:space:]]' ',' | sed 's/,$//g')
CB_FASTQ=$(ls */${SAMPLE_ID}_*R1_001.fastq.gz | tr '[[:space:]]' ',' | sed 's/,$//g')

echo "Aligning samples w/ STARsolo for: ${SAMPLE}."
/home/bnphan/src/STAR-2.7.9a/bin/Linux_x86_64/STAR \
--outFileNamePrefix ${SAMPLE_ID}. \
--soloType Droplet \
--readFilesIn $cDNA_FASTQ $CB_FASTQ \
--readFilesCommand zcat \
--genomeDir $GENOME_DIR \
--limitOutSJcollapsed 5000000 \
--runThreadN 8 \
--soloCBwhitelist $BARCODES \
--soloFeatures GeneFull Gene Velocyto \
--soloBarcodeReadLength 0 \
--soloCBmatchWLtype 1MM \
--soloCellFilter EmptyDrops_CR \
--soloMultiMappers EM \
--soloUMIlen 12 \
--soloUMIdedup 1MM_CR \
--outSAMtype None
# --outSAMtype BAM SortedByCoordinate \
# --outSJtype None \
# --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM CY UY
rsync --remove-source-files -Pauv $TMPDIR/${SAMPLE_ID}* $PROJDIR/data/raw_data/STARsolo_out
rm -rf ${SAMPLE_ID}._STARtmp */${SAMPLE_ID}* ${SAMPLE_ID}*
else 
	echo "A completed STARsolo already exists: ${SOLODIR}/${SAMPLE_ID}.Solo.out"
fi
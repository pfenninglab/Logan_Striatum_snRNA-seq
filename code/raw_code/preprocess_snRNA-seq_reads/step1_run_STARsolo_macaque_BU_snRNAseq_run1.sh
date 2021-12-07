#!/bin/bash
#SBATCH -n 1
#SBATCH --partition=pfen1
#SBATCH --job-name=STARsolo2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=60G
#SBATCH --error=logs/STARsolo2_%A_%a.txt
#SBATCH --output=logs/STARsolo2_%A_%a.txt
#SBATCH --array=1-4

source ~/.bashrc

###########################################
# get sample name in demultiplex csv table
PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=${PROJDIR}/data/raw_data
FASTQDIR=${DATADIR}/BU_snRNAseq_Run1/2021_10_08_LoganR_RayM-295941646/2021_10_08_LoganR_RayM-470955485
SOLODIR=$DATADIR/STARsolo_out

mkdir -p $TMPDIR $SOLODIR; cd $SOLODIR

# sample-specific files
SAMPLE_ID=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $1}' ${PROJDIR}/data/tidy_data/tables/OUD1_snRNA_seq_sampleSheet.csv )
TISSUE=$(awk -F ',' -v IND=${SLURM_ARRAY_TASK_ID} 'FNR == (IND + 1) {print $2}' ${PROJDIR}/data/tidy_data/tables/OUD1_snRNA_seq_sampleSheet.csv )
PREFIX=${SAMPLE_ID}.${TISSUE}.

##################################################
# cDNA fragment in Read2, cell barcode in Read1
cDNA_FASTQ=$(ls ${FASTQDIR}/${SAMPLE_ID}/*R2_001.fastq.gz | tr '[[:space:]]' ',' | sed 's/,$//g')
CB_FASTQ=$(ls ${FASTQDIR}/${SAMPLE_ID}/*R1_001.fastq.gz | tr '[[:space:]]' ',' | sed 's/,$//g')

# perform mapping & read quantification
GENOME_DIR=/home/bnphan/resources/genomes/GRCh38.p13
BARCODES=/home/bnphan/src/cellranger-6.0.0/lib/python/cellranger/barcodes/3M-february-2018.txt

if [[ ! -f "${PREFIX}Log.final.out" ]]; then
	echo "Aligning samples w/ STARsolo for: ${SAMPLE}."
	~/src/STAR-2.7.9a/bin/Linux_x86_64/STAR \
	--outFileNamePrefix $PREFIX \
	--soloType Droplet \
	--readFilesIn $cDNA_FASTQ $CB_FASTQ \
	--readFilesCommand zcat \
	--genomeDir $GENOME_DIR \
	--limitOutSJcollapsed 5000000 \
	--runThreadN 12 \
	--soloCBwhitelist $BARCODES \
	--soloFeatures GeneFull \
	--soloBarcodeReadLength 0 \
	--soloCBmatchWLtype 1MM \
	--soloCellFilter EmptyDrops_CR \
	--soloMultiMappers EM \
	--soloUMIlen 12 \
	--soloUMIdedup 1MM_CR \
	--outSAMtype BAM SortedByCoordinate \
	--outSJtype None \
	--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM CY UY
	rm -rf ${PREFIX}_STARtmp
else 
	echo "A completed STARsolo already exists: ${SOLODIR}/${PREFIX}Solo.out"
fi
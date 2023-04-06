#################
## directories ##
SETDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
GWASDIR=/projects/pfenninggroup/machineLearningForComputationalBiology/gwasEnrichments
CODEDIR=${SETDIR}/code/raw_code/ldsc_interferon_response
DATADIR=${SETDIR}/data/tidy_data/ldsc_interferon_response
ANNOTDIR=${DATADIR}/annotations
cd $CODEDIR; mkdir -p logs $ANNOTDIR; mv -n annot* logs

#######################################
## merge human epigenome backgrounds ##
BG1=${DATADIR}/bed
BG2=/home/bnphan/src/atac-seq-pipeline/genome/hg38/ataqc/reg2map_honeybadger2_dnase_all_p10_ucsc.hg19_to_hg38.bed.gz
BGNAME=IFN_Macrophages_Monocytes.BG_Honeybadger
BGFILE=${DATADIR}/${BGNAME}.bed.gz


checkFile(){
HASFILE=TRUE
for CHR in {1..22}; do F1=$(echo $FILE2 |sed "s/@/$CHR/g");
if [[ ! -f $F1 ]]; then echo "no annotation found for chr${CHR}." && HASFILE=FALSE;
elif [[ $F1 -ot $BED ]]; then
echo "newer foreground file than annotations for chr${CHR}." && HASFILE=FALSE;
fi
done
}

#################################
## make the background peak file
if [[ ! -f ${BGFILE} ]]; then
cat ${BG1}/*.bed.gz | zcat | \
cut -f 1-3 | sort --parallel=10 -k1,1 -k2,2n | gzip > ${DATADIR}/${BGNAME}.tmp.bed.gz 
bedtools merge -i ${DATADIR}/${BGNAME}.tmp.bed.gz | gzip > ${BGFILE}
rm ${DATADIR}/${BGNAME}.tmp.bed.gz
fi


#######################################
# for background for binary annotations
if [[ ! -f "${ANNOTDIR}/${BGNAME}.AFR.1.l2.ldscore.gz" ]]; then 
	sbatch --mem 10G -p pfen2 --job-name=${BGNAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
	-i ${BGFILE} -n ${BGNAME} -g hg38 -p AFR -o $ANNOTDIR
fi
if [[ ! -f "${ANNOTDIR}/${BGNAME}.EUR.1.l2.ldscore.gz" ]]; then 
	sbatch --mem 2G -p pfen2 --job-name=${BGNAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
	-i ${BGFILE} -n ${BGNAME} -g hg38 -p EUR -o $ANNOTDIR
fi


#############################################
## annotate for LDSC w/ binary annotations ##
CTS_AFR_FN=${DATADIR}/IFN_Macrophages_Monocytes.AFR_hg38.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${DATADIR}/IFN_Macrophages_Monocytes.EUR_hg38.ldcts; > $CTS_EUR_FN
for BED in ${BG1}/*.bed.gz ; do
NAME=$(basename $BED .bed.gz)

## check that AFR annotations are complete
FILE2=${ANNOTDIR}/${NAME}.AFR.@.l2.M; checkFile
if [[ ! -f "${ANNOTDIR}/${NAME}.AFR.1.l2.ldscore.gz" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 30G -p pool1 --array=1-22%4 --job-name=${NAME}.AFR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p AFR -o $ANNOTDIR
sleep 10m
fi

## check that AFR annotations are complete
FILE2=${ANNOTDIR}/${NAME}.EUR.@.l2.M; checkFile
if [[ ! -f "${ANNOTDIR}/${NAME}.EUR.1.l2.ldscore.gz" ]]; then 
echo "Annotations for ${NAME} not found."
sbatch --mem 30G -p pool1 --array=1-22%4 --job-name=${NAME}.EUR ${GWASDIR}/scripts/annotate_bed_LDSC_1000G_hg19_hg38.sh \
-i ${BED} -n ${NAME} -g hg38 -p EUR -o $ANNOTDIR
sleep 10m
fi
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.AFR.,${ANNOTDIR}/${BGNAME}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${NAME}.EUR.,${ANNOTDIR}/${BGNAME}.EUR." >> $CTS_EUR_FN
done





##################################################################################
## create an LDCTS file where the IFN has the shared Naive peaks as the background
CTS_AFR_FN=${DATADIR}/IFN_Macrophages_Monocytes_Naive_BG.AFR_hg38.ldcts; > $CTS_AFR_FN
CTS_EUR_FN=${DATADIR}/IFN_Macrophages_Monocytes_Naive_BG.EUR_hg38.ldcts; > $CTS_EUR_FN

for IND in {1..7}; do ## there are 7 datasets w/ treatment and a matching naive bg
TX=$(awk -F"," -v IND=$IND 'NR==(IND + 1){print $7}' $SETDIR/data/tidy_data/ldsc_interferon_response/tables/interferon_macrophages_monocyte_annotations.csv)
NAIVE=$(awk -F"," -v IND=$IND 'NR==(IND + 1){print $8}' $SETDIR/data/tidy_data/ldsc_interferon_response/tables/interferon_macrophages_monocyte_annotations.csv)
NAME=$(basename $TX _Unique)

if [[ $TX != '' ]]; then
echo -e "${NAME}\t${ANNOTDIR}/${TX}.AFR.,${ANNOTDIR}/${NAIVE}.AFR." >> $CTS_AFR_FN
echo -e "${NAME}\t${ANNOTDIR}/${TX}.EUR.,${ANNOTDIR}/${NAIVE}.EUR." >> $CTS_EUR_FN
fi
done




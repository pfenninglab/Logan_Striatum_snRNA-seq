PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
CODEDIR=$PROJDIR/code/final_code/compare_to_Glud1tg_microarray
DATADIR=$PROJDIR/data/tidy_data/compare_to_Glud1tg_microarray

mkdir -p $CODEDIR $DATADIR/sra

cd $DATADIR/sra
tar -xvf GSE11419_RAW.tar # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11419
tar -xvf GSE48911_RAW.tar # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48911

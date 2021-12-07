PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/HeKleyman2021_macaque_striatum_data_processing
mkdir -p $DATADIR/rdas

cd $DATADIR/rdas
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167920/suppl/GSE167920%5FResults%5FMSNs%5Fprocessed%5Ffinal%2Erds%2Egz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE167nnn/GSE167920/suppl/GSE167920%5FResults%5Ffull%5Fnuclei%5Fprocessed%5Ffinal%2Erds%2Egz

gunzip GSE167920_Results_MSNs_processed_final.rds.gz
gunzip GSE167920_Results_full_nuclei_processed_final.rds.gz

PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/Seney2021_OUDvsCTL_bulk_DLPFC_NAc_tables
mkdir -p $DATADIR/tables


cd $DATADIR/tables

wget -O Seney2021_OUD_RNA-seq_DEGs.Data_File_S2.DLFPC.xlsx https://ars.els-cdn.com/content/image/1-s2.0-S000632232101369X-mmc4.xlsx
wget -O Seney2021_OUD_RNA-seq_DEGs.Data_File_S3.NAc.xlsx https://ars.els-cdn.com/content/image/1-s2.0-S000632232101369X-mmc5.xlsx

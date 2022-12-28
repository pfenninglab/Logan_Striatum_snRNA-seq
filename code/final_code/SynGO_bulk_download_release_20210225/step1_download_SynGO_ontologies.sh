PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/SynGO_bulk_download_release_20210225
mkdir -p $DATADIR/tables


cd $DATADIR/tables

wget -O SynGO_bulk_download_release_20210225.zip https://www.syngoportal.org/data/SynGO_bulk_download_release_20210225.zip
unzip SynGO_bulk_download_release_20210225.zip
rm SynGO_bulk_download_release_20210225.zip

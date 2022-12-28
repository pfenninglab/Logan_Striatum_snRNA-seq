PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
DATADIR=$PROJDIR/data/tidy_data/Piechota_et_al_mouse_striatum_on_drugs
mkdir -p $DATADIR/tables

cd $DATADIR/tables
wget --no-check-certificate https://static-content.springer.com/esm/art%3A10.1186%2Fgb-2010-11-5-r48/MediaObjects/13059_2010_2341_MOESM2_ESM.XLS



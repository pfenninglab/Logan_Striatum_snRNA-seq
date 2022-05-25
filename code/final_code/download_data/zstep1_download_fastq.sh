#!/bin/bash 

PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
cd $PROJDIR

mkdir ${PROJDIR}/data/raw_data/BU_snRNAseq_Pilot ${PROJDIR}/data/raw_data/BU_snRNAseq_Run1 ${PROJDIR}/data/raw_data/BU_snRNAseq_Run2 

### download BCL files
bs download run -i 215462301 -o ${PROJDIR}/data/raw_data/BU_snRNAseq_Pilot
bs download run -i 232949763 -o ${PROJDIR}/data/raw_data/BU_snRNAseq_Run1
bs download run -i 233627432 -o ${PROJDIR}/data/raw_data/BU_snRNAseq_Run2
bs download run -i 234243031 -o ${PROJDIR}/data/raw_data/BU_snRNAseq_Run3

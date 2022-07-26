#!/bin/bash

###########################################
######## transfer from Lane to SCC ########
LANEDIR=bnphan@lanec1.compbio.cs.cmu.edu:/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
BU_DIR=/restricted/projectnb/singlecell-rl/Logan_BU_Striatum_snRNA-seq

mkdir -p $BU_DIR/data/raw_data


rsync -Paun $LANEDIR/code $BU_DIR

rsync -Pau $LANEDIR/code $BU_DIR

rsync -Paun $LANEDIR/data/tidy_data $BU_DIR/data

rsync -Pau $LANEDIR/data/tidy_data $BU_DIR/data


###########################################
######## transfer from SCC to Lane ########
rsync -Paun $BU_DIR/* $LANEDIR

rsync -Pau $BU_DIR/* $LANEDIR


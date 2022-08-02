#!/bin/bash

###########################################
######## transfer from Lane to SCC ########
LANEDIR=bnphan@lanec1.compbio.cs.cmu.edu:/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
BU_DIR=/restricted/projectnb/singlecell-rl/Logan_BU_Striatum_snRNA-seq

mkdir -p $BU_DIR/data/raw_data

rsync -Pa $LANEDIR/code $BU_DIR

rsync -Pa $LANEDIR/figures $BU_DIR

rsync -Pa $LANEDIR/data/tidy_data $BU_DIR/data

###########################################
######## transfer from SCC to Lane ########
rsync -Pan $BU_DIR/* $LANEDIR

rsync -Pa $BU_DIR/* $LANEDIR


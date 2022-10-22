#!/bin/bash

###########################################
######## transfer from Lane to SCC ########
LANEDIR=bnphan@lanec1.compbio.cs.cmu.edu:/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
BU_DIR=/restricted/projectnb/singlecell-rl/Logan_BU_Striatum_snRNA-seq

mkdir -p $BU_DIR/data/raw_data

rsync -Pavh --delete $LANEDIR/code $BU_DIR

rsync -Pavh --delete $LANEDIR/figures $BU_DIR

rsync -Pavh --delete $LANEDIR/data/tidy_data $BU_DIR/data

###########################################
######## transfer from SCC to Lane ########
rsync -Pavh $BU_DIR/ $LANEDIR


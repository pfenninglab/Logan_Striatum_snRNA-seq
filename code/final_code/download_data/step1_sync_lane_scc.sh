#!/bin/bash

###########################################
######## transfer from Lane to SCC ########
LANEDIR=bnphan@lanec1.compbio.cs.cmu.edu:/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
BU_DIR=/restricted/project/singlecell-rl/Logan_BU_Striatum_snRNA-seq

mkdir -p $BU_DIR/data/raw_data

# code
rsync -Pav $LANEDIR/code $BU_DIR

# figures
rsync -Pav $LANEDIR/figures $BU_DIR

# products
rsync -Pav $LANEDIR/products $BU_DIR

# tidy_data
rsync -Pav $LANEDIR/data/tidy_data $BU_DIR/data

# raw_data
rsync -Pav $LANEDIR/data/raw_data/STARsolo_out $BU_DIR/data/raw_data


###########################################
######## transfer from SCC to Lane ########
rsync -Paun $BU_DIR/data $LANEDIR
rsync -Paun $BU_DIR/code $LANEDIR
rsync -Paun $BU_DIR/figures $LANEDIR

rsync -Pau $BU_DIR/data $LANEDIR
rsync -Pau $BU_DIR/code $LANEDIR
rsync -Pau $BU_DIR/figures $LANEDIR


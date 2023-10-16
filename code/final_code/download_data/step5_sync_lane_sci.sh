#!/bin/bash

###########################################
######## transfer from Lane to SCI ########
LANEDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
SCI_DIR=badoi.phan3-umw@hpc.umassmed.edu:/pi/ryan.logan-umw/badoi


LANEDIR1=bnphan@lanec1.compbio.cs.cmu.edu:/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq
LANEDIR2=bnphan@lanec1.compbio.cs.cmu.edu:/projects/pfenninggroup/singleCell/Logan_BU_DLPFC_snRNA-seq
LANEDIR3=bnphan@lanec1.compbio.cs.cmu.edu:/projects/pfenninggroup/singleCell/McLean_chronic_opioid_monkey_snRNA-seq
SCI_DIR=/pi/ryan.logan-umw/badoi

rsync -Pavh $LANEDIR1 $SCI_DIR
rsync -Pavh $LANEDIR2 $SCI_DIR
rsync -Pavh $LANEDIR3 $SCI_DIR

###########################################
######## transfer from SCI to Lane ########
rsync -Pavh $BU_DIR/ $LANEDIR


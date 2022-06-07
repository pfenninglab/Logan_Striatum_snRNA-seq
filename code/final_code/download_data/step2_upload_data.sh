#!/bin/bash 

## upload to BU_Addiction_snRNA-seq/Striatum
gdrive mkdir Research/Data/BU_Addiction_snRNA-seq/Striatum_data/STARsolo_out
gdrive sync upload /projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/raw_data/STARsolo_out 177VP-flTiKfYwdA8ZkuijkXoI7cVpjT1

gdrive mkdir Research/Data/BU_Addiction_snRNA-seq/Striatum_data/fastq_merged
gdrive sync upload /projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/raw_data/fastq_merged 1-y6CyfD8yIvpNILW8ibXw83cK-5DGDek

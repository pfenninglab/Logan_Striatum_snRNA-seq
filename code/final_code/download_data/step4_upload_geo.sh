ln -s /data/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/raw_data/fastq_merged \
/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/raw_data/fastq

PROJDIR=/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq

## calculate MD5 Sum
find /data/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/raw_data/fastq_merged -type f | \
parallel -j 8 md5sum > ${PROJDIR}/data/tidy_data/tables/fastq_md5.txt

ls ${PROJDIR}/data/tidy_data/Seurat_projects/BU_OUD_Striatum_refined_all_SeuratObj_N22.* | \
parallel -j 2 md5sum > ${PROJDIR}/data/tidy_data/tables/processed_md5.txt

## send the data to GEO, using user login and password
lftp -u geoftp,'Jirgen5KnyFliefU' ftp://geoftp:********@ftp-private.ncbi.nlm.nih.gov
cd uploads/badoi.phan@pitt.edu_01AX6gqO
mirror -R /data/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/raw_data/fastq_merged
mirror -R /projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/tidy_data/geo_objects/



rclone sync -v /projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/tidy_data/geo_objects/BU_OUD_Striatum_refined_all_SeuratObj_N22.cellxgene.h5ad \
pfenninglab_gdrive:Data/Logan_OUD_Striatum_snRNA/data/tidy_data/geo_objects/

md5sum /projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/tidy_data/geo_objects/BU_OUD_Striatum_refined_all_SeuratObj_N22.cellxgene.h5ad > \
/projects/pfenninggroup/singleCell/Logan_BU_Striatum_snRNA-seq/data/tidy_data/geo_objects/BU_OUD_Striatum_refined_all_SeuratObj_N22.cellxgene.h5ad.md5.txt


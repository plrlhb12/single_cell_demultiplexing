#!/bin/bash/
set -e
set -u
set -o pipefail

# make sure follow the data file structure
# current/sample/
# current/cellsnp.sh
# for the co-runing cellsnp and vireo, use un-subsetted donor file

folder=$1
echo process the sample of $folder
cd $folder
# run cellSNP
[ ! -d "cellSNP_out" ] && mkdir cellSNP_out/ && echo "created cellSNP_out folder"
BAM=possorted_genome_bam.bam
BARCODE=barcodes.tsv.gz
CELLSNP_OUT_DIR=cellSNP_out/
REGION_VCF=../genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
cellSNP -s $BAM -b $BARCODE -O $CELLSNP_OUT_DIR -R $REGION_VCF -p 20 --minMAF 0.1 --minCOUNT 20

# run vireoSNP
[ ! -d "vireo_out" ] && mkdir vireo_out/ && echo "created vireo_out folder"
CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
DONOR_FILE=../merged.vcf.gz
VIREO_OUT_DIR=vireo_out/
D_NUM=10
vireo -c $CELL_FILE -d $DONOR_FILE -N $D_NUM -o $VIREO_OUT_DIR
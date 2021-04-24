#!/bin/bash/

# this is the one run in biowulf

cd cellSNP/
# run cellSNP
mkdir ./cellSNP_out/
BAM=possorted_genome_bam.bam
BARCODE=barcodes.tsv.gz
CELLSNP_OUT_DIR=cellSNP_out/
REGION_VCF=genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz
cellSNP -s $BAM -b $BARCODE -O $CELLSNP_OUT_DIR -R $REGION_VCF -p 20 --minMAF 0.1 --minCOUNT 20

# run vireoSNP
mkdir ./vireo_out/
CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
DONOR_FILE=UNHS_to_use_7s.vcf.gz
VIREO_OUT_DIR=vireo_out/
vireo -c $CELL_FILE -d $DONOR_FILE -o $VIREO_OUT_DIR


#!/bin/bash
set -e
set -u

# run vireoSNP
folder=$1
echo processing the sample of $folder
cd $folder
[ ! -d "vireo_out" ] && mkdir vireo_out/ && echo "created the folder of vireo_out"
CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
DONOR_FILE=sub.merged.vcf.gz
VIREO_OUT_DIR=vireo_out/
D_NUM=10
vireo -c $CELL_FILE -d $DONOR_FILE -N $D_NUM -o $VIREO_OUT_DIR
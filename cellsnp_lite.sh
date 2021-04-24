#!/bin/bash

mkdir cellsnp_lite_out_all_2
time cellsnp-lite -S file_list  -i sample_list -O cellsnp_lite_out_all_2 -p 22 --minMAF 0.1 --minCOUNT 100 --cellTAG None --UMItag None --genotype --gzip
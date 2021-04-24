
# prepare
###################
## copy bam and bai files of the three pools into HPC, notice the missing or "r" in the r3 bai file
## download the individual 10 vcf files from google drive
https://drive.google.com/drive/folders/1dKsBHCMmKucq4XCORh066grW1yLAT7EZ

NIH11A-NIH21A (not include 15A)

## merge individual files into one vcf using bcftools merge, NIH15A is not one of them
## note that in the origal individual single vcf, there isn't chr and pos info. After merge by id, these info show up.

ls *.vcf > file_list

## for h19 vcf
for file in `cat file_list`
do
echo $file
bcftools view --threads 8 $file -Oz -o ${file}.gz;
bcftools index sorted_${file}.gz
done


ls *gz > zip_file_list

## since no chrom and pos infor in vcf, use -m id
bcftools merge --threads 8 -m id -l zip_file_list -Oz -o merged.vcf.gz
bcftools index merged.vcf.gz


# update on hg38_vcf 

for file in `cat file_list`
do
echo $file
bcftools sort $file -Oz -o sorted_${file}.gz;
bcftools index sorted_${file}.gz
done

## since there are chr and pos info in _h38_vcf, skip -m id

bcftools merge --threads 8 -l zip_file_list -Oz -o merged.vcf.gz
bcftools index merged.vcf.gz

## transfer the two files to hpc
## go ahead to prepare cellsnp

## prepare barcodes.tsv.gz from the .h5 file
## notice that the adata.obs is still an empty df, the index is actually the values of barcodes

write a script to retrieve barcodes and zip 
get_barcodes.py
```
import scanpy as sc
adata = sc.read_10x_h5()
adata.obs.to_csv('NIHPlR1.tsv', sep = '\t', header=None)

gzip NIHPlR1.tsv

```
## transfer the NIHPlR1.tsv.gz, NIHPlR2.tsv.gz, NIHPlR3.tsv.gz files to hpc
scp -r NIH*.tsv.gz pengl7@helix.nih.gov:/data/CARD/tprojects/projects/Xylena/

###### make the file structure as sample/cellSNP, sample/cellsnp.sh
mkdir -p r{1,2,3}

#### tee ~/folder1/test.txt ~/folder2/test.txt < ~/test.txt

for num in 1 2 3
do
folder1=NIHPlR${num}
folder2=NIHPlR${num}/cellSNP
mkdir -p ${folder2}
cp /data/CARD/tprojects/projects/cellSNP/Combined2_forebrain/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz ${folder2}/
cp /data/CARD/tprojects/projects/cellSNP/Combined2_forebrain/cellsnp.sh ${folder1}
cp merged.vcf.gz* ${folder2}/
mv ${folder1}possorted_genome_bam.bam* ${folder2}/
mv ${folder1}.tsv.gz ${folder2}/
done

for num in 1 2 3
do
folder1=NIHPlR${num}
folder2=NIHPlR${num}/cellSNP
rm ${folder1}/cellsnp.sh
echo "removing old cellsnp.sh"
cp /data/CARD/tprojects/projects/cellSNP/Combined2_forebrain/cellsnp_2021.sh ${folder1}/cellsnp_2021.sh
echo "copy new cellsnp.sh"
mv ${folder2}/${folder1}possorted_genome_bam.bam ${folder2}/possorted_genome_bam.bam
mv ${folder2}/${folder1}possorted_genome_bam.bam.bai ${folder2}/possorted_genome_bam.bam.bai
mv ${folder2}/${folder1}.tsv.gz ${folder2}/barcodes.tsv.gz
done

for num in 1 2 3
do
folder1=NIHPlR${num}
folder2=NIHPlR${num}/cellSNP
mv ${folder1}/cellsnp.sh ${folder1}/cellsnp_2021.sh
done

<!-- ## modify cellsnp.swarm file
## failed becasue timeout
swarm -f cellsnp_2021.swarm -g 120 -t 48 --time 24:00:00 --module cellsnp/0.1.7,vireosnp/0.3.2

6386557
sjobs
sacct

## rerun 
jobs were time out after running 24hrs.
The useages of cpu and mem show only use 16g mem and 21 cpu 
find out that inside the cellsnp_2021.sh file, -p 20; change to 40, however forgot to change it iside each sample folder. !!!

swarm -f cellsnp_2021.swarm -g 20 -t 40 --time 48:00:00 --module cellsnp/0.1.7,vireosnp/0.3.2

6461934

the run of r2 was killed after cellsnp
The memory use exceeded memory allocation, and one or more processes were killed

re-run the vireo part and give more cpu and mem
can define the number of donors -N 10

swarm -f cellsnp_2021.swarm -g 60 -t 40 --time 24:00:00 --module cellsnp/0.1.7,vireosnp/0.3.2


## received the data from pool3 again, redo anaylsis, using -p 40

swarm -f cellsnp_2021.swarm -g 20 -t 40 --time 48:00:00 --module cellsnp/0.1.7,vireosnp/0.3.2

6524082

but stuck at the last step of cellsnp forever -->

#######################
#### 21Jan 9:00pm again

`swarm -f cellsnp_2021.swarm -g 360 -t 40 --time 56:00:00 --module cellsnp/0.1.7,vireosnp/0.3.2`

6609127

# remember to add cd cellSNP/ in the vireo_2021.sh file

<!-- The memory use exceeded memory allocation, and one or more processes were killed

also killed because format not the same

Traceback (most recent call last):
  File "/opt/conda/envs/app/bin/vireo", line 8, in <module>
    sys.exit(main())
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/vireo.py", line 142, in main
    cell_dat, donor_vcf = match_donor_VCF(cell_dat, donor_vcf)
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py", line 20, in match_donor_VCF
    if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:
TypeError: '>=' not supported between instances of 'NoneType' and 'int' -->

##  
for whatever reasons : the generated cellSNP.cells.vcf.gz has the "chr" prefix while the merged donor and genom1k all just have numeric chrosomes

my previous sequ all doesn't have chr

`bcftools query -f '%CHROM\n' NIHPlR1/cellSNP/cellSNP_out/cellSNP.cells.vcf.gz`

chr1
chr2

this format is decided by the chr format in bam file: 
`samtools view NIHPlR1/possorted_genome_bam.bam | grep 'chr' | less`

A00428.... chr1 ....


`bcftools query -f '%CHROM\n' NIHPlR1/cellSNP/merged.vcf.gz`

1
2
3

`bcftools query -f '%CHROM\n' NIHPlR1/cellSNP/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz`

1
2
3

so change chr on merged vcf files

`cp ../compare-variants/chr_name_conv.txt ./`

`bcftools annotate --threads 56 --rename-chrs chr_name_conv.txt  merged.vcf.gz -Oz -o merged_chr.vcf.gz`
`bcftools index --threads 56 merged_chr.vcf.gz`



`$ vireo -c $CELL_FILE -d $DONOR_FILE -o $VIREO_OUT_DIR`

<!-- [pengl7@cn3300 cellSNP]$ vireo -c $CELL_FILE -d $DONOR_FILE -o $VIREO_OUT_DIR
[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py:17: RuntimeWarning: invalid value encountered in greater_equal
  if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:
/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py:20: RuntimeWarning: invalid value encountered in greater_equal
  if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:
[vireo] 87 out 224764 variants matched to donor VCF
[vireo] Demultiplex 14249 cells to 1 donors with 87 variants.
[vireo] lower bound ranges [-8932.3, -8932.3, -8932.3]
[vireo] allelic rate mean and concentrations:
[[0.271 0.276 0.28 ]]
[[5366.7 5366.7 5366.7]]
[vireo] donor size before removing doublets:
donor0
14249
Traceback (most recent call last):
  File "/opt/conda/envs/app/bin/vireo", line 8, in <module>
    sys.exit(main())
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/vireo.py", line 187, in main
    ASE_mode=options.ASE_mode)
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/vireo_wrap.py", line 122, in vireo_wrap
    doublet_prob, ID_prob = modelCA.predict_doublet(AD, DP)
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/vireo_model.py", line 313, in predict_doublet
    GT_both = add_doublet_GT(self.GT_prob)
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/vireo_model.py", line 394, in add_doublet_GT
    s_idx1 = sp_idx[:, 0]
IndexError: too many indices for array -->

#####################
## use conda env to run vireo separately
`sinteractive --mem=240g --cpus-per-task=20`

cd cellSNP/

## made some subset on the donor file

wc -l merged.chr.vcf.gz
19249 merged.chr.vcf.gz

wc -l NIHPlR1/cellSNP/cellSNP_out/cellSNP.cells.vcf.gz
250034 NIHPlR1/cellSNP/cellSNP_out/cellSNP.cells.vcf.gz

`bcftools view merged.chr.vcf.gz -R NIHPlR1/cellSNP/cellSNP_out/cellSNP.cells.vcf.gz -Oz -o NIHPlR1/sub.merged.vcf.gz`

less NIHPlR1/sub.merged.vcf.gz
less NIHPlR2/sub.merged.vcf.gz
less NIHPlR3/sub.merged.vcf.gz

for folder in NIHPlR3
do
bcftools view merged.chr.vcf.gz -R ${folder}/cellSNP_out/cellSNP.cells.vcf.gz -Oz -o ${folder}/sub.merged.vcf.gz
done

wc -l r1.sub.merged.vcf.gz
423 r1.sub.merged.vcf.gz

wc -l NIHPlR1/cellSNP/sub.merged.vcf.gz
16 NIHPlR1/cellSNP/sub.merged.vcf.gz

CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
IREO_OUT_DIR=vireo_out/
DONOR_FILE=sub.merged.vcf.gz

vireo -c $CELL_FILE -d $DONOR_FILE -N 11 -o $VIREO_OUT_DIR

1:30pm

<!-- [pengl7@cn3300 cellSNP]$ vireo -c $CELL_FILE -d $DONOR_FILE -N 11 -o $VIREO_OUT_DIR
[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py:17: RuntimeWarning: invalid value encountered in greater_equal
  if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:
/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py:20: RuntimeWarning: invalid value encountered in greater_equal
  if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:
[vireo] 87 out 224764 variants matched to donor VCF
[vireo] Demultiplex 14249 cells to 11 donors with 87 variants.
[vireo] lower bound ranges [-6483.2, -6236.0, -5476.1]
[vireo] allelic rate mean and concentrations:
[[0.033 0.283 0.881]]
[[5552.2 8423.5 2124.3]]
[vireo] donor size before removing doublets:
donor0	donor1	donor2	donor3	donor4	donor5	donor6	donor7	donor8	donor9	donor10
1283	1290	1298	1293	1293	1291	1297	1311	1305	1292	1297
[vireo] final donor size:
unassigned
14249
Traceback (most recent call last):
  File "/opt/conda/envs/app/bin/vireo", line 8, in <module>
    sys.exit(main())
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/vireo.py", line 201, in main
    donor_GPb[idx, :, :], donor_vcf['samples'])
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/plot/base_plot.py", line 44, in plot_GT
    import matplotlib.pyplot as plt
  File "/opt/conda/envs/app/lib/python3.7/site-packages/matplotlib/__init__.py", line 207, in <module>
    _check_versions()
  File "/opt/conda/envs/app/lib/python3.7/site-packages/matplotlib/__init__.py", line 192, in _check_versions
    from . import ft2font
ImportError: /opt/conda/envs/app/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so: failed to map segment from shared object -->

failed again: but this time has output, takes about 40 mins without plotting

## rerun without giving donor file and and number of donors
vireo -c $CELL_FILE -o $VIREO_OUT_DIR
#not working 

sinteractive --mem=240g --cpus-per-task=20
source /data/$USER/conda/etc/profile.d/conda.sh
conda env list
conda activate 

#############################
## run without donor file but defining number of donors
## saved in 22Jan, vireo_out_noDonor
vireo -c $CELL_FILE -N 10 -o $VIREO_OUT_DIR

1hr50imn + 1 hr50min not finish the step of plotting: abut 4 hours, 
successfuly generated 10 donors, but failed to match them to the 10 source donors

##############################
25Jan2021
# swarm to finish the cellsnp step for NIHPlR3 sample
`swarm -f cellsnp_2021_up.swarm -g 40 -t 20 --time 48:00:00 --module cellsnp/0.1.7`

6695817

run 1d + 19 hours for NIHPlR3 to finish the cellsnp step


######### wait for the new vcf files after liftOver
######### 

`swarm -f vireo_2021_up.swarm -g 360 -t 20 --time 8:00:00 --module vireosnp/0.3.2`

<!-- 7246853
waiting time is too long due to big mem, cancel the job -->

swarm -f vireo_2021_up.swarm --partition='quick' -g 360 -t 20 --time 3:59:00 --module vireosnp/0.3.2

7269882

# success partially three runs in this batch started at 30Jan, 12:00
# all used subsettted donor files and version 0.3.2
7279138

it took 1 hour to generate all asssign but failed at plotting due to errror in im

[+] Loading vireosnp  0.3.2  on cn2402 
[+] Loading singularity  3.7.1  on cn2402 
/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py:17: RuntimeWarning: invalid value encountered in greater_equal
  if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:
/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py:20: RuntimeWarning: invalid value encountered in greater_equal
  if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:
Traceback (most recent call last):
  File "/opt/conda/envs/app/bin/vireo", line 8, in <module>
    sys.exit(main())
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/vireo.py", line 201, in main
    donor_GPb[idx, :, :], donor_vcf['samples'])
  File "/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/plot/base_plot.py", line 44, in plot_GT
    import matplotlib.pyplot as plt
  File "/opt/conda/envs/app/lib/python3.7/site-packages/matplotlib/__init__.py", line 207, in <module>
    _check_versions()
  File "/opt/conda/envs/app/lib/python3.7/site-packages/matplotlib/__init__.py", line 192, in _check_versions
    from . import ft2font
ImportError: /opt/conda/envs/app/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so: failed to map segment from shared object
swarm_7279138_1.e (END)

check the donor_id.txt, quite weired results, about 5 thousands unassigned due to co-assign and also others.

### 
/usr/local/lmod/modulefiles

using version 0.4.2
test by adjustting vireo_run_up.sh using unsubsetted merage filee and not using D-NUM in vireo_run_up_2.sh

fisrt mv out output

for folder in NIHPlR1 NIHPlR2 NIHPlR3
do
mv ${folder}/vireo_out ${folder}/vireo_out_1
done


/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP

30Jan
# only try NIHPlR3 first, pending for several hours so cancel it , should specificy using quick partition
freen # even quick has no resource

<!-- # use NIHPlR2 and norm module to test

# using unsubseted donor file
swarm -f vireo_2021_up.swarm -g 240 -t 20 --time 2:00:00 --module vireosnp/0.3.2

'''
CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
DONOR_FILE=./merged.chr.vcf.gz
VIREO_OUT_DIR=vireo_out/
D_NUM=11
vireo -c $CELL_FILE -d $DONOR_FILE -N D_NUM -o $VIREO_OUT_DIR
'''
didn't find the output folder
7308565 -->

# at the samee time usying my conda install new version to test NIHPlR1
## must specify the number of donors, otherwise two many indexes for array


import cellsnp

sinteractive --mem=240g --cpus-per-task=20
source /data/pengl7/conda/etc/profile.d/conda.sh
conda activate vireo-env
cd NIHPlR1

## 1. test1 witthout specify -N
mkdir vireo_0.4_out/

CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
VIREO_OUT_DIR=vireo_0.4_out/
DONOR_FILE=../merged.chr.vcf.gz

2:06pm
2:49
vireo -c $CELL_FILE -d $DONOR_FILE -o $VIREO_OUT_DIR

[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
[vireo] 5179 out 224764 variants matched to donor VCF
[vireo] Demultiplex 14249 cells to 1 donors with 5179 variants.
[vireo] lower bound ranges [-520996.1, -520996.1, -520996.1]
[vireo] allelic rate mean and concentrations:
[[0.353 0.353 0.353]]
[[343082.7 343082.7 343082.7]]
[vireo] donor size before removing doublets:
donor0
14249

Traceback (most recent call last):
  File "/data/pengl7/conda/envs/vireo-env/bin/vireo", line 8, in <module>
    sys.exit(main())
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/vireo.py", line 187, in main
    ASE_mode=options.ASE_mode)
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/utils/vireo_wrap.py", line 122, in vireo_wrap
    doublet_prob, ID_prob = modelCA.predict_doublet(AD, DP)
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/utils/vireo_model.py", line 313, in predict_doublet
    GT_both = add_doublet_GT(self.GT_prob)
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/utils/vireo_model.py", line 394, in add_doublet_GT
    s_idx1 = sp_idx[:, 0]
IndexError: too many indices for array: array is 1-dimensional, but 2 were indexed

## test2. -N 10 , donor=subsetted
mkdir vireo_0.4_out_n10_sub/

CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
VIREO_OUT_DIR=vireo_0.4_out_n10_sub/  # orignially vireo_0.4_out_2/
DONOR_FILE=sub.merged.vcf.gz

vireo -c $CELL_FILE -d $DONOR_FILE -N 10 -o $VIREO_OUT_DIR

[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
[vireo] 5179 out 224764 variants matched to donor VCF
[vireo] Demultiplex 14249 cells to 10 donors with 5179 variants.
[vireo] lower bound ranges [-308607.9, -268972.4, -260671.0]
[vireo] allelic rate mean and concentrations:
[[0.045 0.455 0.894]]
[[401111.  493929.6 134207.4]]
[vireo] donor size before removing doublets:
donor0	donor1	donor2	donor3	donor4	donor5	donor6	donor7	donor8	donor9
600	2575	2088	1144	627	1299	1672	1448	903	1893
[vireo] final donor size:
NIH_11A	NIH_12A	NIH_13A	NIH_14A	NIH_16A	NIH_17A	NIH_18A	NIH_19A	NIH_20A	NIH_21A	doublet	unassigned
156	1847	1417	701	308	817	1062	601	283	975	654	5428
Traceback (most recent call last):
  File "/data/pengl7/conda/envs/vireo-env/bin/vireo", line 8, in <module>
    sys.exit(main())
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/vireo.py", line 201, in main
    donor_GPb[idx, :, :], donor_vcf['samples'])
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/plot/base_plot.py", line 45, in plot_GT
    import matplotlib.pyplot as plt
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/__init__.py", line 174, in <module>
    _check_versions()
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/__init__.py", line 159, in _check_versions
    from . import ft2font
ImportError: /data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so: failed to map segment from shared object: Cannot allocate memory

## test3, n=11, donor=unsubsetted
mkdir vireo_0.4_out_n11/
CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
VIREO_OUT_DIR=vireo_0.4_out_n11/  # origianl vireo_0.4_out/
DONOR_FILE=../merged.chr.vcf.gz

vireo -c $CELL_FILE -d $DONOR_FILE -N 11 -o $VIREO_OUT_DIR

less than 1 hour

(vireo-env) [pengl7@cn3284 Xylena]$ cd NIHPlR1
(vireo-env) [pengl7@cn3284 NIHPlR1]$ ELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
(vireo-env) [pengl7@cn3284 NIHPlR1]$ CELL_FILE=cellSNP_out/cellSNP.cells.vcf.gz
(vireo-env) [pengl7@cn3284 NIHPlR1]$ VIREO_OUT_DIR=vireo_0.4_out_n11/
(vireo-env) [pengl7@cn3284 NIHPlR1]$ DONOR_FILE=../merged.chr.vcf.gz
(vireo-env) [pengl7@cn3284 NIHPlR1]$ vireo -c $CELL_FILE -d $DONOR_FILE -N 11 -o $VIREO_OUT_DIR
[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
[vireo] 5179 out 224764 variants matched to donor VCF
[vireo] Demultiplex 14249 cells to 11 donors with 5179 variants.
[vireo] lower bound ranges [-308650.6, -265778.3, -257220.7]
[vireo] allelic rate mean and concentrations:
[[0.04  0.452 0.899]]
[[390476.9 507473.4 131297.7]]
[vireo] donor size before removing doublets:
donor0	donor1	donor2	donor3	donor4	donor5	donor6	donor7	donor8	donor9	donor10
511	688	1983	2488	1209	1096	675	635	1555	606	2804
[vireo] final donor size:
NIH_11A	NIH_12A	NIH_13A	NIH_14A	NIH_16A	NIH_17A	NIH_18A	NIH_19A	NIH_20A	NIH_21A	donor1	doublet	unassigned
65	209	1389	1823	784	695	191	192	1035	310	1974	597	4985
Traceback (most recent call last):
  File "/data/pengl7/conda/envs/vireo-env/bin/vireo", line 8, in <module>
    sys.exit(main())
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/vireo.py", line 201, in main
    donor_GPb[idx, :, :], donor_vcf['samples'])
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/plot/base_plot.py", line 45, in plot_GT
    import matplotlib.pyplot as plt
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/__init__.py", line 174, in <module>
    _check_versions()
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/__init__.py", line 159, in _check_versions
    from . import ft2font
ImportError: /data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so: failed to map segment from shared object: Cannot allocate memory


########## # 1Feb2021
# 1. after Jean Mao install vireo 0.4.0 , this is done on NIHPlR2 on Feb 1, 5:43pm
swarm -f vireo_2021_up.swarm -g 240 -t 20 --time 2:00:00 --partition 'quick' --module vireosnp/0.4.2
7444880
still error

`swarm -f vireo_2021_up.swarm -g 240 -t 20 --time 2:00:00 --partition 'quick' --module vireosnp/0.4.2`

# 2. test old file from iNDI
`swarm -f vireo.swarm -g 240 -t 20 --time 2:00:00 --partition 'quick' --module vireosnp/0.3.2`
7454302
no error at tall

# 3. test wehter the chr in the orignal cellSNP.cells.vcf.gz caused trouble:

processing the sample of Combined2_forebrain
created the folder of vireo_out
[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
[vireo] 41391 out 69192 variants matched to donor VCF
[vireo] Demultiplex 9064 cells to 10 donors with 41391 variants.
[vireo] lower bound ranges [-1378268.3, -1270974.4, -1224304.1]
[vireo] allelic rate mean and concentrations:
[[0.02  0.472 0.969]]
[[1675585.  1998158.3  839761.7]]
[vireo] donor size before removing doublets:
donor0  donor1  donor2  donor3  donor4  donor5  donor6  donor7  donor8  donor9
759     497     1288    1936    925     268     2085    286     344     676
[vireo] final donor size:
GT19-38445      GT19-38446      GT19-38447      GT19-38448      GT19-38449      GT19-38450      GT19-38451      donor7  donor8  donor9  doublet unassigned
687     454     1149    1811    852     249     1939    136     291     523     695     278
[vireo] All done: 26 min 21.9 sec


although still have error message as 
[+] Loading vireosnp  0.3.2  on cn2278
[+] Loading singularity  3.7.1  on cn2278
/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py:17: RuntimeWarning: invalid value encountered in greater_equal
  if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:
/opt/conda/envs/app/lib/python3.7/site-packages/vireoSNP/utils/io_utils.py:20: RuntimeWarning: invalid value encountered in greater_equal
  if np.sum(mm_idx == mm_idx) == 0 or np.sum(mm_idx >= 0) == 0:


# turns out not

change chr1 to 1

echo "1 chr1" >> chr_name_conv.txt
echo "2 chr2" >> chr_name_conv.txt


cp NIHPlR1/cellSNP_out/ NIHPlR1/cellSNP_out_old
cp chr_name_conv.txt NIHPlR1/cellSNP_out/

change column seq
awk '{ print $2 " " $1}' chr_name_conv.txt chr_name_conv_res.txt

cann't use bcftools --rename becasue chr is not defineed in header
bcftools annotate --rename-chrs chr_name_conv_res.txt cellSNP.cells.vcf.gz | bgzip > renamed.cellSNP.cells.vcf.gz


awk '{gsub(/^chr/,""); print}' cellSNP.base.vcf > no_chr.cellSNP.base.vcf

for file in cellSNP.cells.vcf
do
gunzip ${file}.gz;
sed 's/chr//g' $file > no_chr.$file;
gzip no_chr.${file};
mv no_chr.${file}.gz ${file}.gz
done

Feb1, 22:50pm
use merged file withou chr to do it again
the cellsnp_out file was saved as cellsnp_out_change_chr
change on vireo_2021_up_3.sh: change merge file to the one using the format of "1"
save the output in vireo_out_3_change_chr (orginal vireo_out then vireo_out_3)

`swarm -f vireo_2021_up_3.swarm -g 240 -t 20 --time 2:00:00 --partition 'quick' --module vireosnp/0.3.2`

7467293

the same error

3Feb2021: solutions
1. using foreceLearnGT in vireo
2. redo-cellsnp using h19 reference
3. re-generate donor vcf using cellsnp-lite, note that 

Note that now it is recommended to use the version of reference with AF5e2 (i.e., AF>5%, 7.4M SNPS), instead of AF5e4 (i.e., AF>0.05%, 36.6M SNPs).

lists of the reference snp datasets
https://sourceforge.net/projects/cellsnp/files/SNPlist/



###### read vireo manual again, try the option of -foreLearnGT

activate the vireo-env

cd ../NIHPlR1
mkdir vireo_0.4_forceLearnGT

CELL_DATA=cellSNP_out/cellSNP.cells.vcf.gz
DONOR_GT_FILE=sub.merged.vcf.gz
OUT_DIR=vireo_0.4_forceLearnGT

vireo -c $CELL_DATA -d $DONOR_GT_FILE -o $OUT_DIR -N 10 --forceLearnGT --randSeed 2

# no much improvement

[vireo] Loading cell VCF file ...

[vireo] Loading donor VCF file ...
[vireo] 5179 out 224764 variants matched to donor VCF
[vireo] Demultiplex 14249 cells to 10 donors with 5179 variants.
[vireo] lower bound ranges [-324019.1, -269557.6, -263628.8]
[vireo] allelic rate mean and concentrations:
[[0.044 0.451 0.893]]
[[393229.9 502445.1 133572.9]]
[vireo] donor size before removing doublets:
donor0	donor1	donor2	donor3	donor4	donor5	donor6	donor7	donor8	donor9
713	1307	714	1170	2568	1905	1312	2066	791	1701
[vireo] final donor size:
NIH_11A	NIH_12A	NIH_13A	NIH_14A	NIH_16A	NIH_17A	NIH_18A	NIH_19A	NIH_20A	NIH_21A	doublet	unassigned
113	817	319	709	1846	1035	553	1404	181	1063	653	5556




## check the number of variants of donor file
wc -l /Users/pengl7/Downloads/Xylena/deconvolution/vcfs/h19/NIH_11A_Variants.vcf 
487412 /Users/pengl7/Downloads/Xylena/deconvolution/vcfs/h19/NIH_11A_Variants.vcf


wc -l merged_copy.vcf
408015 merged_copy.vcf

wc -l sub.merged.vcf
8042 sub.merged.vcf

wc -l UNHS_to_use_7s.vcf
90870 UNHS_to_use_7s.vcf

41391 out 69192 variants matched to donor VCF

wc -l UNHS_to_use.vcf 
   90868 UNHS_to_use.vcf

grep 'indel' merged_copy.vcf
grep 'indel' merged_copy.vcf | wc -l
10512

$ grep 'indel' /Users/pengl7/Downloads/Xylena/deconvolution/vcfs/h19/NIH_11A_Variants.vcf | wc -l
11028


408015
 10512


# check homozyous 
bcftools view -i 'GT[*]="alt"' merged_copy.vcf > heter_merged.vcf

wc -l heter_merged.vcf
  208455 heter_merged.vcf


bcftools view -i 'GT[*]="alt"' UNHS_to_use.vcf > homo_unhs.vcf

 wc -l homo_unhs.vcf
   90870 homo_unhs.vcf

# not work for merged vcf, becasue no label of snp
bcftools view heter_merged.vcf -i 'TYPE="snp"' > snp_heter_merged.vcf

# suggestions from Yuanhua
Hi Lirong,

If you have separate scRNA for each individual, then using these scRNA for genotyping would give you good overlapped SNPs. You can use cellSNP to perform genotyping in a bulk manner, for each individual separately or jointly (mode 2b):

https://cellsnp-lite.readthedocs.io/en/latest/manual.html#mode-2-pileup-whole-chromosome-s-without-given-snps

Note, this is re-implemented cellSNP with C/C++, achieving ~10 times faster.
please add --cellTAG None --UMItag None --genotype, but don't provide cell barcode -b as you treat them as a bulk sample. You can use prefiltering, e.g., --minMAF 0.1 --minCOUNT 100.

Yuanhua


NIH12Apossorted_genome_bam.bam

$OUT_DIR
sinteractive --mem=40g --cpus-per-task=20

mkdir cellsnp_lite_out

time cellsnp-lite -s NIH12Apossorted_genome_bam.bam -I NIH12A -O cellsnp_lite_out -p 22 --minMAF 0.1 --minCOUNT 100 --cellTAG None --UMItag None --genotype --gzip


[I::main] All Done!
[I::main] end time: 2021-02-05 10:46:20
[I::main] time spent: 2419 seconds.

real	40m19.637s
user	437m54.999s
sys	0m56.293s
ls *bam > file_list
sed 's/possorted_genome_bam.bam//g' file_list > sample_list

mkdir cellsnp_lite_out_all
time cellsnp-lite -S file_list  -i sample_list -O cellsnp_lite_out_all -p 22 --minMAF 0.1 --minCOUNT 100 --cellTAG None --UMItag None --genotype --gzip

Killed due to mem over 40 and cpu over 20

real	177m12.401s
user	3485m51.727s
sys	10m37.738s


# generate cellsnp_lite.sh and cellsnp_lite.swarm

swarm -f cellsnp_lite.swarm -g 240 -t 40 --time 24:00:00 --module cellsnp-lite

7746910

finished in 12hr10min

[I::main] All Done!
[I::main] end time: 2021-02-06 03:20:33
[I::main] time spent: 43815 seconds.

real    730m22.373s
user    6578m4.393s
sys     39m46.626s

[pengl7@biowulf individuals]$ wc -l cellsnp_lite_out_all_2/cellSNP.cells.vcf.gz
790854 cellsnp_lite_out_all_2/cellSNP.cells.vcf.gz

redo vireo now
change the name of donor file
cp cellsnp_lite_out_all_2/cellSNP.cells.vcf.gz ../cellSNP.donor.vcf.gz


note that that this donor file's chrom using numberical, so have to change it. However, bcftools annotate doesn't work because %CHROM isn't in the header.

so change the cellSNP.cells.vcf.gz

for folder in NIHPlR{1,2,3}
do
#gunzip ${folder}/cellSNP_out/cellSNP.cells.vcf.gz
#sed 's/chr//g' ${folder}/cellSNP_out/cellSNP.cells.vcf > ${folder}/cellSNP_out/nochr_cellSNP.cells.vcf
bgzip ${folder}/cellSNP_out/cellSNP.cells.vcf
#bgzip ${folder}/cellSNP_out/nochr_cellSNP.cells.vcf
done



##

change vireo_2021_up_3.sh

run NIHRPl1 and 2 first using swarm

swarm -f vireo_2021_up_3.swarm -t 240 -g 20 --time 3:00:00 --module vireosnp
no free node

test NIHPRl1 in conda env
source /data/pengl7/conda/etc/profile.d/conda.sh

mkdir vireo_out_6Feb2021
CELL_FILE=cellSNP_out/nochr_cellSNP.cells.vcf.gz
DONOR_FILE=../cellSNP.donor.vcf.gz
OUT_DIR=vireo_out_6Feb2021
time vireo -c $CELL_FILE -d $DONOR_FILE -N 10 -o $OUT_DIR --noPlot

1:00pm


[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
[vireo] 128186 out 224764 variants matched to donor VCF
[vireo] Demultiplex 14249 cells to 10 donors with 128186 variants.
[vireo] lower bound ranges [-6007259.3, -6007259.3, -6007259.3]
[vireo] allelic rate mean and concentrations:
[[0.054 0.452 0.876]]
[[7301663.3 7243071.3 1637751.4]]
[vireo] donor size before removing doublets:
donor0	donor1	donor2	donor3	donor4	donor5	donor6	donor7	donor8	donor9
1190	140	1302	541	1720	402	3060	3247	2344	303
[vireo] final donor size:
NIH11A	NIH12A	NIH13A	NIH14A	NIH16A	NIH17A	NIH18A	NIH19A	NIH20A	NIH21A	doublet	unassigned
852	107	996	318	1322	284	2188	2509	1777	176	2595	1125
Traceback (most recent call last):
  File "/data/pengl7/conda/envs/vireo-env/bin/vireo", line 8, in <module>
    sys.exit(main())
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/vireo.py", line 203, in main
    plot_GT(out_dir, res_vireo['GT_prob'][idx, :, :], donor_names)
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/vireoSNP/plot/base_plot.py", line 45, in plot_GT
    import matplotlib.pyplot as plt
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/__init__.py", line 174, in <module>
    _check_versions()
  File "/data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/__init__.py", line 159, in _check_versions
    from . import ft2font
ImportError: /data/pengl7/conda/envs/vireo-env/lib/python3.7/site-packages/matplotlib/ft2font.cpython-37m-x86_64-linux-gnu.so: failed to map segment from shared object: Cannot allocate memory

real	42m59.958s
user	34m46.585s
sys	6m56.958s
(vireo-env) [pengl7@cn3167 NIHPlR1]$



(vireo-env) [pengl7@cn3167 NIHPlR2]$ mkdir vireo_out_6Feb2021
(vireo-env) [pengl7@cn3167 NIHPlR2]$ CELL_FILE=cellSNP_out/nochr_cellSNP.cells.vcf.gz
(vireo-env) [pengl7@cn3167 NIHPlR2]$ DONOR_FILE=../cellSNP.donor.vcf.gz
(vireo-env) [pengl7@cn3167 NIHPlR2]$ OUT_DIR=vireo_out_6Feb2021
(vireo-env) [pengl7@cn3167 NIHPlR2]$ time vireo -c $CELL_FILE -d $DONOR_FILE -N 10 -o $OUT_DIR --noPlot
[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
[vireo] 128681 out 227449 variants matched to donor VCF
[vireo] Demultiplex 14893 cells to 10 donors with 128681 variants.
[vireo] lower bound ranges [-7074553.4, -7074553.4, -7074553.4]
[vireo] allelic rate mean and concentrations:
[[0.066 0.447 0.849]]
[[8427696.9 8307120.4 1862431.8]]
[vireo] donor size before removing doublets:
donor0	donor1	donor2	donor3	donor4	donor5	donor6	donor7	donor8	donor9
1204	186	1270	568	1813	611	3037	3470	2444	291
[vireo] final donor size:
NIH11A	NIH12A	NIH13A	NIH14A	NIH16A	NIH17A	NIH18A	NIH19A	NIH20A	NIH21A	doublet	unassigned
816	125	874	293	1303	439	2105	2532	1688	155	2831	1732
[vireo] All done: 45 min 23.2 sec


real	46m5.364s
user	37m2.934s
sys	7m35.004s

for NIHPRl3 it tooks 320G and more than 2hrs

swarm -f vireo_2021_up.swarm -g 370 -t 10 --time 3:59:00 --partition 'quick' --module vireosnp

7812951

finished in 55 mins
max mem use 326g, 2cpu

processing the sample of NIHPlR3
created the folder of vireo_out_out
[vireo] Loading cell VCF file ...
[vireo] Loading donor VCF file ...
[vireo] 153748 out 274302 variants matched to donor VCF
[vireo] Demultiplex 16867 cells to 10 donors with 153748 variants.
[vireo] lower bound ranges [-9121672.0, -9121672.0, -9121672.0]
[vireo] allelic rate mean and concentrations:
[[0.069 0.444 0.844]]
[[10635556.8 10621582.9  2384148.3]]
[vireo] donor size before removing doublets:
donor0  donor1  donor2  donor3  donor4  donor5  donor6  donor7  donor8  donor9
1384    195     1370    654     1971    471     3808    3899    2780    334
[vireo] final donor size:
NIH11A  NIH12A  NIH13A  NIH14A  NIH16A  NIH17A  NIH18A  NIH19A  NIH20A  NIH21A  doublet unassigned
990     151     997     403     1474    329     2603    2885    2099    200     3312    1424



### do some cleaning, for example, remove bam files, remove intermediate files, change intermediate file name to conversion file name, delete useless results

for folder in NIHPlR{1,2,3}
do
ls -lth ${folder}/*Feb*
done


for folder in NIHPlR{1,2,3}
do
ls -lth ${folder}
done

for folder in NIHPlR{1,2,3}
do
mkdir -p demultiplexing_out_new/${folder}
cp -r ${folder}/cellSNP_out demultiplexing_out_new/${folder}/cellSNP_out
cp -r ${folder}/vireo_out_6Feb2021 demultiplexing_out_new/${folder}/vireo_out_6Feb2021
done

for folder in NIHPlR{1,2,3}
do
ls -lth demultiplexing_out_new/${folder}/vireo_out_6Feb2021
ls -lth demultiplexing_out_new/${folder}/cellSNP_out
done


for folder in NIHPlR{1,2,3}
do
rm demultiplexing_out_new/${folder}/cellSNP_out/cellSNP.cells.vcf
mv demultiplexing_out_new/${folder}/cellSNP_out/nochr_cellSNP.cells.vcf.gz demultiplexing_out_new/${folder}/cellSNP_out/cellSNP.cells.vcf.gz
done


for folder in NIHPlR{1,2,3}
do
rm demultiplexing_out_new/${folder}/cellSNP_out/cellSNP.cells.vcf
mv demultiplexing_out_new/${folder}/cellSNP_out/nochr_cellSNP.cells.vcf.gz demultiplexing_out_new/${folder}/cellSNP_out/cellSNP.cells.vcf.gz
done


# change back the file name and delete intermediate files

for folder in NIHPlR{1,2,3}
do
rm ${folder}/cellSNP_out/cellSNP.cells.vcf
mv ${folder}/cellSNP_out/nochr_cellSNP.cells.vcf.gz ${folder}/cellSNP_out/cellSNP.cells.vcf.gz
done


for folder in NIHPlR{1,2,3}
do
echo $folder >> summary_all3.tsv
cat ${folder}/vireo_out_6Feb2021/summary.tsv >> summary_all3.tsv
echo  >> summary_all3.tsv
done


NIHPlR1
Var1    Freq
NIH11A  852
NIH12A  107
NIH13A  996
NIH14A  318
NIH16A  1322
NIH17A  284
NIH18A  2188
NIH19A  2509
NIH20A  1777
NIH21A  176
doublet 2595
unassigned      1125

NIHPlR2
Var1    Freq
NIH11A  816
NIH12A  125
NIH13A  874
NIH14A  293
NIH16A  1303
NIH17A  439
NIH18A  2105
NIH19A  2532
NIH20A  1688
NIH21A  155
doublet 2831
unassigned      1732

NIHPlR3
Var1    Freq
NIH11A  990
NIH12A  151
NIH13A  997
NIH14A  403
NIH16A  1474
NIH17A  329
NIH18A  2603
NIH19A  2885
NIH20A  2099
NIH21A  200
doublet 3312
unassigned      1424

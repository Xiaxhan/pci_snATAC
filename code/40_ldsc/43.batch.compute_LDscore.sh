#!/bin/bash
#conda activate ldsc
indir=outMarkerP_LDscore

path_ldsc=ldsc
path_data=ldsc/data
path_resource=ldsc_ph/resources

group=$1
chr=$2

echo ">$group.$chr"
             bfile=$path_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr
             outname=$group.$chr
             snpfile=hm.$chr.snp
             python $path_ldsc/ldsc.py \
                     --l2 \
                     --bfile $bfile \
                     --ld-wind-cm 1 \
                     --annot $path_data/$indir/$outname.annot.gz \
                     --thin-annot \
                     --out $path_data/$indir/$outname \
                     --print-snps $path_resource/hapmap3_snps/$snpfile

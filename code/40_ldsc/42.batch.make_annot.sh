#!/bin/bash
#conda activate ldsc

outdir=outMarkerP_LDscore
indir=dxMarkerPeak

path_ldsc=ldsc
path_resource=ldsc_ph/resources
path_data=ldsc/data

path_outAnno=$path_data/$outdir

chrs=`seq 1 22`

group=$1
chr=$2
	bed=$path_data/$indir/$group"_dxMarkerP_ctUniqueP.bed"
	ls $bed
		annot_file=$group.$chr.annot.gz
		bim=$path_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr.bim
		echo $annot_file
		python $path_ldsc/make_annot.py --bed-file $bed --bimfile $bim --annot-file $path_outAnno/$annot_file 




#!/bin/bash
#conda activate allcools
chrom_size=sorted.hg38.chrom.sizes
mC_path=03_bgzip_tbi

bed_path=Dx_peaks_bed

input_mC=sort.hg38.allc_Astro.tsv.gz
oprefix=${input_mC#*.hg38.}
oprefix=${oprefix%.tsv*}
echo $oprefix

declare -A bed_paths
declare -A group_names
id=0
for i in `cat list.tbi`
do
        path=$bed_path"/$i"
        name=${i%.sorted*}
        bed_paths[$id]=$path
        group_names[$id]=$name
        id=`expr $id + 1`
done
echo ${bed_paths[*]}
echo ${group_names[*]}

allcools allc-to-region-count --allc_path $mC_path/$input_mC\
        --chrom_size_path $chrom_size\
        --output_prefix $oprefix\
        --mc_contexts CNN\
	--save_zero_cov\
	--cov_cutoff 9999\
	--region_bed_paths ${bed_paths[*]}\
	--region_bed_names ${group_names[*]}

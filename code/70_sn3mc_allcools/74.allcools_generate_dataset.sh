#!/bin/bash
#conda activate allcools
chrom_size=sorted.hg38.chrom.sizes

bed_path=Dx_peaks_bed

bed_gz=bvFTD.odc.C7.sorted.bed.gz
bedname=peak
echo $bedname

allcools generate-dataset\
	--allc_table test_allc_table.tsv\
	--output_path /geschwindlabshares/RexachGroup/Xia_Data/heterchromatin/scripts/DxPeak_mC.mcds\
	--chrom_size_path $chrom_size\
	--obs_dim cell \
	--cpu 2\
	--chunk_size 1\
	--regions $bedname $bed_path/$bed_gz\
	--quantifiers $bedname hyper-score CNN,CHN cutoff=0.1\
	--quantifiers $bedname hypo-score CNN,CHN cutoff=0.1\
	--quantifiers $bedname count CNN,CHN 

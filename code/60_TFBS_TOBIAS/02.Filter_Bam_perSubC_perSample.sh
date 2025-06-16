#!/bin/bash
subc=$1
inp="fragments"

for sample in `cat list.86`
do
	in_barcode="/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/bam/barcodes/${subc}.${sample}.barcodes.txt" 
	SAM_body="/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/bam/sam/body_${sample}.sam" 
	SAM_header="/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/bam/sam/header_${sample}.sam"
	
	# output files
	tmp_clean_sam="/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/bam/${subc}/tmp_${subc}_${sample}.sam"
	filtered_sam="/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/bam/${subc}/${subc}_${sample}.filtered.sam"
	filtered_bam="/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/bam/${subc}/${subc}_${sample}.filtered.bam"
	
	echo "> Extracting reads of $subc from sample $sample"
	
	if [ -f $in_barcode ]; then	
		# Filter alignment and keep subC barcodes in each sample
		echo "--extracting barcodes"
		LC_ALL=C grep -F -f $in_barcode $SAM_body > $tmp_clean_sam

		# Combine header and body
		echo "--combining header and subc sam"
		cat $SAM_header $tmp_clean_sam > $filtered_sam
	
		# Convert filtered.sam to BAM format
		echo "--sam to bam"
		samtools view -b $filtered_sam > $filtered_bam
	
        	echo "Finish BAM: $filtered_bam"
	fi
done

#/bin/bash
outdir="sam"

inp="fragments"

for sample in `cat list.86`
do
	bam="$inp/$sample/outs/possorted_bam.bam"
	if [ ! -f $outdir/header_${sample}.sam ]; then
		echo "Processing sample: $sample"
        	echo "Input BAM: $bam"
    		samtools view -H "$bam" > "$outdir/header_${sample}.sam"   # Extract header
		samtools view "$bam" > "$outdir/body_${sample}.sam"  
	        echo "--Done"	
	fi	

done

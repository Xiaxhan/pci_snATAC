#ATAC correct Tobias - subsample and then perform tn5 shift
#conda activate tobias_env

INDIR=bam
GENOME=GRCh38.primary_assembly.genome.fa
OUTDIR=subc_specific_peaks

PEAKDIR=ldsc/data/cluster_specific_atacPeak

groups=("ast.C1" "mg.C4")
for group in ${groups[@]}
do
	prefix=${group//./} # astC1 mgC4
	bam="Subtype_"$prefix"_sorted.bam"	
    	correct_bam=$OUTDIR/"Subtype_$group"."_corrected.bw"
    	PEAK=$PEAKDIR/$group"_newRaw.bed"
	echo $group

	if [[ ! -f $correct_bam ]]; then
    	TOBIAS ATACorrect --bam $INDIR/$bam \
		--genome $GENOME --peaks $PEAK \
		--outdir $OUTDIR --cores 8
	fi
    
done



#ATAC correct Tobias - subsample and then perform tn5 shift
#conda activate tobias_env

OUTDIR=/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/tobias/subc_specific_peaks

groups=("ast.C1" "mg.C4")
for group in ${groups[@]}
do
	PEAK=/geschwindlabshares/RexachGroup/Xia_Data/ldsc/data/cluster_specific_atacPeak/$group"_newRaw.bed"
	
	prefix=${group//./} # astC1 mgC4
	in_bw="$OUTDIR/Subtype_$prefix""_sorted_corrected.bw"
	out_bw="$OUTDIR/$prefix.footprint.bw"
   	
	if [[ ! -f $out_bw ]]; then
    	TOBIAS FootprintScores --signal $in_bw \
    	--regions $PEAK \
    	--output $out_bw --cores 4
    	fi
done

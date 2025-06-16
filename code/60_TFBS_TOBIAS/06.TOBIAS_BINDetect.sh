#ATAC correct Tobias - subsample and then perform tn5 shift
#conda activate tobias_env

OUTDIR=subc_specific_peaks
GENOME=GRCh38.primary_assembly.genome.fa
MOTIF=motif_databases/CIS-BP_2.00/Homo_sapiens.meme
groups=("ast.C1" "mg.C4")
for group in ${groups[@]}
do
	prefix=${group//./} # astC1 mgC4
	PEAK=/geschwindlabshares/RexachGroup/Xia_Data/ldsc/data/cluster_specific_atacPeak/$group"_newRaw.bed"

	TOBIAS BINDetect --motifs $MOTIF \
                 --signals $OUTDIR'/'$prefix.footprint.bw \
                 --genome $GENOME \
                 --peaks $PEAK \
		 --cond_names $prefix \
                 --outdir $OUTDIR/$prefix --cores 4 \
                 --prefix $prefix
done

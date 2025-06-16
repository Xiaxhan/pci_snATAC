#!/bin/bash
#conda activate ldsc

indir=outMarkerP_LDscore
outdir=outMarkerP_h2


path_ldsc=ldsc
path_data=ldsc/data
path_resource=ldsc_ph/resources
path_gwas=resources/GWAS
BASELINE_1000G=$path_resource/1000G_EUR_Phase3_baseline/baseline.
WEIGHTS_BASE=$path_resource/weights_hm3_no_hla/weights.
FREQ_BASE=$path_resource/1000G_Phase3_frq/1000G.EUR.QC.

declare -A sumstats
sumstats=(
[ADHD_Demontis_ng2023]=ADHD_Demontis_ng2023.sumstats.gz
[ASD_Grove_ng2019]=ASD_Grove_ng2019.sumstats.gz
[MDD_Howard_nn2019]=MDD_Howard_nn2019.sumstats.gz
[SCZ_Trubetskoy_n2022]=SCZ_Trubetskoy_n2022.sumstats.gz
[MS_Andlauer_SciAdv2016]=MS_Andlauer_SciAdv2016.sumstats.gz
[ALS_Rheenen_ng2021]=ALS_Rheenen_ng2021.sumstats.gz
[IBD_Glanville_2021]=IBD_Glanville_2021.sumstats.gz
[AD_Bellenguez_ng2022]=AD_Bellenguez_ng2022.sumstats.gz
[AD_Bellenguez]=AD_Bellenguez_ng2022.sumstats.gz
[AD_Jansen]=AD_Jansen_ng2019.sumstats.gz
[AD_Kuncle]=AD_Kuncle_ng2019.sumstats.gz
[AD_Lambert]=AD_Lambert_ng2013.sumstats.gz
[AD_MorenoGrau]=AD_MorenoGrau_2021.sumstats.gz
[PSP_Hoglinger]=PSP_Hoglinger_ng2011.sumstats.gz
[FTD_bvFTD]=FTD_bvFTD_Ferrari_2014.sumstats.gz
[FTD_MMD]=FTD_MMD.sumstats.gz
[FTD_META]=FTD_META.sumstats.gz
[FTD_PNFA]=FTD_PNFA.sumstats.gz
[FTD_SD]=FTD_SD.sumstats.gz       
)

group=$1
dx=$2

echo ">$group|$dx"

		out=$group.$dx 
		
		file_sumstat=$path_gwas/${sumstats[$dx]}
		python $path_ldsc/ldsc.py \
			--h2 $file_sumstat \
			--out $path_data/$outdir/$out \
			--ref-ld-chr $path_data/$indir/$group.,$BASELINE_1000G  \
			--w-ld-chr $WEIGHTS_BASE \
			--overlap-annot \
			--frqfile-chr $FREQ_BASE \
			--print-coefficients

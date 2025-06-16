#!/bin/bash
chr=$1

outdir=QTL_fsVars/out/mpra_LDblock
bfile=resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr

snp_file="mpra.$chr.txt"
out_prefix="LD_R2_mpra.$chr"

plink --bfile $bfile \
        --r2 \
        --ld-snp-list $outdir/$snp_file \
        --ld-window-kb 1000 \
        --ld-window 99999 \
        --ld-window-r2 0.2 \
        --out $outdir/$out_prefix

#!/bin/bash
chr=$1

bfile=resources/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr
outdir=QTL_fsVars/out/mpra_LDblock

snp_file="mpra.$chr.txt"
out_prefix="LDblock_mpra.$chr"

plink --bfile $bfile \
	--blocks no-pheno-req \
	--extract $outdir/$snp_file \
	--out $outdir/$out_prefix

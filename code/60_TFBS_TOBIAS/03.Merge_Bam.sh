#!/bin/bash

subc=$1
outdir="bam"
merged_bam="Subtype_$subc.bam"


samtools merge $outdir/$merged_bam $outdir/${subc}/*filtered.bam
samtools sort -o "$outdir/${merged_bam%.bam}_sorted.bam" $outdir/$merged_bam
samtools index "$outdir/${merged_bam%.bam}_sorted.bam"

echo "Final merged BAM: $outdir/${merged_bam%.bam}_sorted.bam"

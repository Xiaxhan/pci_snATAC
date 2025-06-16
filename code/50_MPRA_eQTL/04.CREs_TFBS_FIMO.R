# Predict the TFBS of each CRE peak
library(rhdf5)
library(data.table)
library(dplyr)
library(parallel)
library(preprocessCore)
library(BSgenome.Hsapiens.UCSC.hg38)
library(memes)

source("pci_snATAC_peakSet_w_Enh.R")
pciATAC_peakset$peakloci <- paste(pciATAC_peakset@seqnames, 
                                  pciATAC_peakset@ranges@start,
                                  pciATAC_peakset@ranges@start + pciATAC_peakset@ranges@width - 1,
                                  sep = "-") 
dt_peakinfo <- values(pciATAC_peakset[,c("peakloci","peakType","nearestGene")])
rownames(dt_peakinfo) <- dt_peakinfo$peakloci

CREs <- dt_peakinfo[dt_peakinfo$peakType %in% c("PE","Promoter"),] %>% unique()

tmpV <- strsplit(CREs$peakloci, split = "-") %>% unlist()
gr_CREs <- GRanges(
  seqnames=tmpV[seq(1, length(tmpV)-2, by = 3)],
  IRanges(start=as.numeric(tmpV[seq(2, length(tmpV)-1, by = 3)]),
          end=as.numeric(tmpV[seq(3, length(tmpV), by = 3)])
  )
)

# annotate TFBS for peaks given TFs - by MEME's FIMO
PATHDB<- "motif_databases"
motifdb <- file.path(PATHDB, "CIS-BP_2.00/Homo_sapiens.meme")
outdir <- "motif_by_memes"

setwd(outdir)

gr_CREs_seq <- get_sequence(gr_CREs,BSgenome.Hsapiens.UCSC.hg38)
writeXStringSet(gr_CREs_seq, filepath = "all_CREs_267322.fa")
fimo_command <- paste(
  paste("fimo --oc fimo_output_allcres",motifdb,"all_CREs_267322.fa")
)

# run the FIMO command
system(fimo_command)

# load 'motif occurence' results
dt_motif <- read.csv(file.path(outdir,"fimo_output_allcres/fimo.tsv"), sep = "\t")
dt_motif <- dt_motif[!is.na(dt_motif$q.value),]
my_motif <- dt_motif[dt_motif$motif_alt_id %in% TFs,]
max(my_motif$p.value)
gr_motif <- GRanges(seqnames = my_motif$sequence_name,
                    ranges = IRanges(as.numeric(my_motif$start),
                                     as.numeric(my_motif$stop)),
                    TF=my_motif$motif_alt_id,
                    motif_match=paste(my_motif$sequence_name,":",
                                      my_motif$start,"-",
                                      my_motif$stop, sep = ""))

#output bed for IGV
export(upDAR_novel, con = file.path(plotdir, bedfile), format = "BED")
for (tf in TFs){
  bedfile <- paste("UpNoveldaEnh_x",tf,".bed", sep = "")
  out_bed <- file.path(plotdir,".bed")
  export(gr_motif[gr_motif$TF %in% tf,],
         con = file.path(plotdir, bedfile), format = "BED")
}

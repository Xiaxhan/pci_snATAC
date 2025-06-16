library(ggplot2)
library(dplyr)
library(data.table)
library(GenomicRanges)
outp <- "QTL_fsVars/out"
majorCs <- c("ast","mg","EX","IN","odc","opc")
Dxs <- c("Control","AD","bvFTD","PSP_S")

##define MPRA fsVar by cutoff: P < 0.05
mpra_in_celltypePeaks <- read.csv(file.path(outp, "mpra_in_celltypePeaks.csv"), row.names = 1)
mpra_in_celltypePeaks <- mpra_in_celltypePeaks[!is.na(mpra_in_celltypePeaks$chr_hg38),]
fsVar <- mpra_in_celltypePeaks[mpra_in_celltypePeaks$pval < 0.05,] 
fsVar_dyn <- fsVar[fsVar$label %in% "Dynamic", ]

nb_p_g <- NULL
for (ct in majorCs){
  dt <- fsVar_dyn[fsVar_dyn$ct %in% ct,]
  peaks <- unique(dt$loci)
  genes <- unique(dt$nearestGene) %>% .[!is.na(.)]
  dte <- dt[dt$peakType %in% "Enhancer",]
  enh <- unique(dte$loci)
  peg <- unique(dte$PE_G) %>% .[!is.na(.)]
  
  nb_p_g <- rbind(nb_p_g,
             data.frame(ct, 
             npeak = length(peaks),
             nNearestGene=length(genes),
             nEnh = length(enh),
             nPEG=length(peg)
             ))
}
nb_p_g

nb_Dx_pg <- NULL
for (dx in Dxs){
  dt <- fsVar_dyn[!is.na(fsVar_dyn[[dx]]),]
  peaks <- unique(dt$loci)
  genes <- unique(dt$nearestGene) %>% .[!is.na(.)]
  dte <- dt[dt$peakType %in% "Enhancer",]
  enh <- unique(dte$loci)
  peg <- unique(dte$PE_G) %>% .[!is.na(.)]

  ng_ct <- NULL
  for (ct in majorCs){
    dct <- dt[dt$PE_G %in% peg & dt$ct %in% ct, ]
    ng_ct <- c(ng_ct, length(unique(dct$PE_G)))
  }
  names(ng_ct) <- majorCs
  row_add <- cbind(data.frame(dx, 
                        npeak = length(peaks),
                        nNearestGene=length(genes),
                        nEnh = length(enh),
                        nPEG=length(peg)), t(as.data.frame(ng_ct)))
  nb_Dx_pg <- rbind(nb_Dx_pg, row_add)
}
nb_Dx_pg

gr <- fsVar_dyn[fsVar_dyn$ct %in% "mg" & fsVar_dyn$peakType %in% "Enhancer",]
gr <- unique(gr$loci)
tmp <- strsplit(gr, split = ":") %>% unlist()
gr <- GRanges(seqnames = tmp[seq(1, length(tmp)-2, by = 3)],
                  ranges = IRanges(as.numeric(tmp[seq(2, length(tmp)-1, by = 3)]), 
                                   end = as.numeric(tmp[seq(3, length(tmp), by = 3)]) 
                                   ))
library(memes)
library(universalmotif)
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
sequences <- get_sequence(gr, genome)

outdir <- "/home/xhan/data_cellrangerOut/project0318/projrmSubset_subPeak/motif"
PATHDB<- "/geschwindlabshares/RexachGroup/SharedData/motif_databases"
motifdb <- file.path(PATHDB, "CIS-BP_2.00/Homo_sapiens.meme")
memeBin <- "/home/xhan/apps/meme-5.4.1/bin"

# motif discovery for given sequences
dir <- file.path(outdir, "Mg_27dynEnh_wMPRAfsVar")
if (!dir.exists(dir)){
  dir.create(dir)
}
findMotif <- runMeme(sequences, 
                     outdir = dir,
                     meme_path = memeBin,
                     parse_genomic_coord = FALSE)
setwd(dir)
system(paste("tomtom -no-ssc -oc . -verbosity 1 -min-overlap 5 -mi 1 -dist pearson -evalue -thresh 10.0 -time 300 meme.txt", motifdb))


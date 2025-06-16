library(biomaRt)
library(data.table)
library(dplyr)
library(GenomicRanges)
p1 <- "resources"
mpra_pos <- read.csv(file.path(p1,"MPRA_w_pos_GRCh38.csv"))

# load atac peaks
path_mypeak <- "QTL_fsVars/mypeaks"
atac_cts <- c("ast","mg","EX","IN","odc","opc")
bed_files <- read.csv(file.path(path_mypeak,"list"), header = F)$V1
list_atacpeaks <- list()
for (file in bed_files){
  fname <- gsub(".bed","",file)
  fname <- gsub("_dxInvariantP_ctUniqueP", ".invariantP",fname)
  fname <- gsub("_dxMarkerP_ctUniqueP", ".dynamicP",fname)
  fname <- gsub(".unique_peaks", "",fname)
  
  dt <- read.table(file.path(path_mypeak, file))
  gr <- GRanges(seqnames = dt$V1, ranges = IRanges(dt$V2, end = dt$V3))
  
  # add PE-gene info
  m <- findMatches(gr, pciATAC_peakset)
  gr$peakType <- rep(NA, length(gr))
  gr$peakGene <-  rep(NA, length(gr))
  gr[queryHits(m)]$peakType <- as.character(pciATAC_peakset[subjectHits(m)]$peakType)
  gr[queryHits(m)]$peakGene <- as.character(pciATAC_peakset[subjectHits(m)]$nearestGene)
  
  list_atacpeaks[[fname]] <- gr
}
names(list_atacpeaks)
names(list_atacpeaks)[names(list_atacpeaks) %like% "EX"]

Get_Peak_wSNPs <- function(dt_snp){
  dt_overlap <- NULL
  gr_qtl <- GRanges(seqnames = paste("chr",dt_snp$chr_name,sep = ""), 
                    ranges = IRanges(dt_snp$chrom_start, end = dt_snp$chrom_start), 
                    snp=dt_snp$refsnp_id)
  # four peak sets: eg "ast.invariantP" "ast.dynamicP"   "ast.Distal"     "ast.PE" ...
  peak_sets <- names(list_atacpeaks)
  for (peakset in peak_sets){
    gr_peak <- list_atacpeaks[[peakset]]
    o <- findOverlaps(gr_peak, gr_qtl)
    
    tmp_peak <- as.data.frame(gr_peak)
    gr_peak$loci <- paste(tmp_peak$seqnames, tmp_peak$start, tmp_peak$end, sep = ":")
    
    dt_merge <- cbind(data.frame(gr_peak[queryHits(o)]), 
                        data.frame(gr_qtl[subjectHits(o)]))
    dt_merge$peakset <- rep(peakset, dim(dt_merge)[1])
    
    df.overlap <- rbind(df.overlap, dt_merge)
  }
  
  return(df.overlap)
}

peaks_w_snps <- Get_Peak_wSNPs(mpra_pos)

write.csv(peaks_w_snps, file.path(outp, "MPRA_snps_in_eachpeak.csv"))
mygenes <- c("STX7","STX4","CSK","NOTCH1","AFG1L",
             "SOX10","PLP1","PLLP")

peaks_w_snps[peaks_w_snps$peakGene %in% mygenes, ]

nb_sum <- NULL
for (peakT in unique(peaks_w_snps$peakset)){
  dt <- peaks_w_snps[peaks_w_snps$peakset %in% peakT,]
  nb_g <- length(unique(dt$peakGene))
  nb_snps <- length(unique(dt$snp))
  dts <- data.frame(peakT, nb_gene=nb_g, nb_snp=nb_snps)
  nb_sum <- rbind(nb_sum, dts)
}
nb_sum

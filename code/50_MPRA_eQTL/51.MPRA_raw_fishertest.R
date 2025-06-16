library(biomaRt)
library(data.table)
library(dplyr)
library(GenomicRanges)
listEnsembl()
outp <- "/geschwindlabshares/RexachGroup/Xia_Data/QTL_fsVars/out"

source("/geschwindlabshares/RexachGroup/Xia_Data/My_R_libs/pci_snATAC_peakSet_w_Enh.R")
colors.ct <- c(ArchR::ArchRPalettes$stallion[c(1,5,3,11)], "#b9a802","#308ac4")
names(colors.ct) <- c("ast","mg","IN","EX","odc","opc")

fetch_snp_pos <- function(snps){
  ensembl_snps <- useEnsembl(biomart="snps", dataset="hsapiens_snp")
  snp_pos <- getBM(attributes = c('refsnp_id','allele','chr_name','chrom_start','chrom_strand', 'minor_allele_freq'), 
                   filters = c('snp_filter'), 
                   values = snps,
                   mart = ensembl_snps)
  snp_pos <- unique(snp_pos)
  snp_pos <- snp_pos[!snp_pos$chr_name %like% "[A-Za-z]", ]
  rownames(snp_pos) <- snp_pos$refsnp_id
  return(snp_pos)
}

count_Peak_wSNPs <- function(dt_snp){
  dt_overlap <- NULL
  gr_qtl <- GRanges(seqnames = paste("chr",dt_snp$chr_name,sep = ""), 
                    ranges = IRanges(dt_snp$chrom_start, end = dt_snp$chrom_start), 
                    snp=dt_snp$refsnp_id)
  peak_sets <- names(list_atacpeaks)
  for (peakset in peak_sets){
    gr_peak <- list_atacpeaks[[peakset]]
    o <- findOverlaps(gr_peak, gr_qtl)
    
    tmp_peak <- as.data.frame(gr_peak)
    gr_peak$loci <- paste(tmp_peak$seqnames, tmp_peak$start, tmp_peak$end, sep = ":")
    
    df.overlap <- data.frame(gr_peak[queryHits(o)], gr_qtl[subjectHits(o)])
    df.overlap <- df.overlap[,c("loci","peakType","peakGene","snp")]
    df.overlap <- df.overlap[!is.na(df.overlap$peakGene), ]
    nb_peak_w_qtl <- length(unique(df.overlap$loci))
    nb_peak_wo_qtl <- length(gr_peak) - nb_peak_w_qtl
    dt <- data.frame(peakT=peakset,
                     QTLtotal=length(gr_qtl),
                     QTLs_in_peak = length(unique(df.overlap$snp)), 
                     Peaks_w_qtl = nb_peak_w_qtl, 
                     Peaks_wo_qtl = nb_peak_wo_qtl)
    dt_overlap <- rbind(dt_overlap, dt)
  }
  
  return(dt_overlap)
}

run_fisherTest <- function(df, peakA, peakB){
  fisher_out <- NULL
  # make contingency table -
  #if row#1 vs row#2 has a OR on col#1 (given col1, row1 has a OR when compared with row2)
  label <- paste(peakA, "vs", peakB, sep = " ")
  df$RAWpeakT <- df$peakT
  for (ct in c("ast","mg","IN","EX","odc","opc")){
    dtss <- df[df$RAWpeakT %like% ct,]
    dtss$peakT <- gsub(paste(ct,".",sep = ""),"", dtss$peakT)
    rownames(dtss) <- dtss$peakT
    dynamicP_w_qtl <- dtss[peakA,"Peaks_w_qtl"]
    dynamicP_wo_qtl <- dtss[peakA,"Peaks_wo_qtl"]
    invarP_w_qtl <- dtss[peakB,"Peaks_w_qtl"]
    invarP_wo_qtl <- dtss[peakB,"Peaks_wo_qtl"]
    contingency_table <- matrix( c(dynamicP_w_qtl, dynamicP_wo_qtl,
                                   invarP_w_qtl,invarP_wo_qtl), nrow=2,
				byrow=TRUE,
                                 dimnames = list(PeakType = c(peakA, peakB),
                                                 QTL = c("With","Without")))
    print(contingency_table)
    fisher_test <- fisher.test(contingency_table)
    
    dtout <- data.frame(ct, label,
                        OR=fisher_test$estimate, 
                        pval=fisher_test$p.value)
    fisher_out <- rbind(fisher_out, dtout)
  }
  rownames(fisher_out) <- NULL
  return(fisher_out)
}


#load mpra
p1 <- "resources"
mpra <- read.csv(file.path(p1, "MPRA_list_of_variants.tsv"), sep = ",", row.names = 1)

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


#----------
# fetch mpra position (hg38)
mpra_pos <- fetch_snp_pos(mpra$variant_id)
peaks_x_mpra <- count_Peak_wSNPs(mpra_pos)

# fisher test
fisher_CRE <- run_fisherTest(peaks_x_mpra, "CRE", "Distal")
fisher_dy <- run_fisherTest(peaks_x_mpra, "dynamicP", "invariantP")
fisher_dyEnh <- run_fisherTest(peaks_x_mpra, "Enhancer.dynamicP", "Enhancer.invariantP")
fisher_dy_dyCRE <- run_fisherTest(peaks_x_mpra, "CRE.dynamicP", "CRE.invariantP")
fisher_dy_woDistal <- run_fisherTest(peaks_x_mpra, "wo_distal.dynamicP", "wo_distal.invariantP")

fisher_out<- rbind(fisher_CRE, fisher_dy,fisher_dyEnh, fisher_dy_dyCRE,fisher_dy_woDistal)
fisher_out

#csv to excel
write.csv(peaks_x_mpra, file.path(outp, "MPRA_in_peaksets.csv"))
write.csv(fisher_out, file.path(outp, "fishertest.MPRA_in_peakset.csv"))

#plot
fisher_out <- read.csv(file.path(outp, "fishertest.MPRA_in_peakset.csv"), row.names = 1)

unique(fisher_out$label)
update_label <- data.frame(
  raw=unique(fisher_out$label)[-1],
  new=c("All Peaks",
        "Enhancer",
        "Enhancer + Promoter",
        "Enhancer + Promoter + Intronic + Exonic")
)
rownames(update_label) <- update_label$raw
fisher_out$type <- update_label[fisher_out$label,]$new
fisher_out <- fisher_out[!is.na(fisher_out$type),]

fisher_out$p_flag <- ifelse(fisher_out$pval < 0.05, "*", NA)

fisher_out$type <- factor(fisher_out$type, 
                         levels = c("All Peaks",
                                    "Enhancer + Promoter + Intronic + Exonic",
                                    "Enhancer + Promoter","Enhancer"))
p <- ggplot(fisher_out,  aes(x = OR, y = -log10(pval))) + 
  geom_point(aes(color = ct, shape=type), size=3) +
  labs(title="MPRA variant Enrichment in Dynamic/Invariant Peaks ", 
       x="Odds Ratio",y="Fisher Exact p-value (-log10)")+
  xlim(c(0.5,2)) +
  geom_vline(xintercept = 1, linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") + 
  scale_color_manual(values = colors.ct) + 
  theme_set(theme_bw())+
  theme(  plot.title = element_text(size=9, face = "bold"),
          axis.title.x=element_text(size=9,face = "plain"),
          title =element_text(size=9, face='bold'),
          axis.title.y=element_text(size=9,vjust = 2,hjust = 0.5,face = "plain"),
          axis.text.x=element_text(size=9, face="plain", angle = 90, hjust = 1, vjust=0.5),
          #axis.text.x=element_blank(),
          axis.text.y=element_text(size=9,face="plain"),
          legend.position = "right",
          legend.title=element_blank(),
          legend.text=element_text(size=9),
          panel.grid = element_blank(),
          #panel.grid.minor.y = element_line(linewidth=1, colour = "grey"),
          #panel.border = element_rect(linetype = "dashed", fill = NA),
          panel.border = element_blank(),
          axis.ticks = element_line(colour = "black", linewidth = 0.1),
          axis.line = element_line(linewidth = 0.1, colour = "black", linetype="solid"))
pdf(file.path(outp, "enrich.MPRA_fisherTest.pdf"), width = 6, height = 4)
p
dev.off()



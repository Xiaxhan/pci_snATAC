# keep QTLs whose eGene match to Enhancer/promoter-linked Gene
library(data.table)
library(dplyr)
library(GenomicRanges)
source("pci_snATAC_peakSet_w_Enh.R")
colors.ct <- c(ArchR::ArchRPalettes$stallion[c(1,5,3,11)], "#b9a802","#308ac4")
names(colors.ct) <- c("ast","mg","IN","EX","odc","opc")

outp <- "/geschwindlabshares/RexachGroup/Xia_Data/QTL_fsVars/out"

path_mypeak <- "/geschwindlabshares/RexachGroup/Xia_Data/QTL_fsVars/mypeaks"
atac_cts <- c("ast","mg","EX","IN","odc","opc")
bed_files <- read.csv(file.path(path_mypeak,"list"), header = F)$V1
# load atac peaks
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


# load eQTL 
eQTL_Malhotra_w_pos <- readRDS("/geschwindlabshares/RexachGroup/Xia_Data/QTL_fsVars/out/eQTL_Malhotra_wpos.rds")

# sn-eQTL
count_Peak_wQTL <- function(dt_QTL, QTLsource){
  QTL_cts <- unique(dt_QTL$celltype)
  dt_overlap <- NULL
  for (ct in QTL_cts){
    atac_ct <- ct_match[ct,]$atac_ct
    if(!is.na(atac_ct)){
      dt_qtl <- dt_QTL[dt_QTL$celltype %in% ct, ]
      gr_qtl <- GRanges(seqnames = dt_qtl$chr, 
                        ranges = IRanges(dt_qtl$start, end = dt_qtl$start), 
                        snp=dt_qtl$snp, QTLgene = dt_qtl$gene)
      # four peak sets: eg "ast.invariantP" "ast.dynamicP"   "ast.Distal"     "ast.PE" ...
      peak_sets <- names(list_atacpeaks)[names(list_atacpeaks) %like% atac_ct]
      for (peakset in peak_sets){
        gr_peak <- list_atacpeaks[[peakset]]
        o <- findOverlaps(gr_peak, gr_qtl)
        
        tmp_peak <- as.data.frame(gr_peak)
        gr_peak$loci <- paste(tmp_peak$seqnames, tmp_peak$start, tmp_peak$end, sep = ":")
     
        df.overlap <- data.frame(gr_peak[queryHits(o)], gr_qtl[subjectHits(o)])
        df.overlap <- df.overlap[,c("loci","peakType","peakGene","snp","QTLgene")]
        df.overlap <- df.overlap[!is.na(df.overlap$peakGene), ]
        #  count overlap only when eGene=peakGene
        keep.overlap <- df.overlap[df.overlap$peakGene == df.overlap$QTLgene,]
        nb_peak_w_qtl <- length(unique(keep.overlap$loci))
        nb_peak_wo_qtl <- length(gr_peak) - nb_peak_w_qtl
        peakT <- gsub(paste(atac_ct,".",sep = ""),"",peakset)
        dt <- data.frame(atac_ct, peakT, QTL_ct=ct, QTLsource, 
                         QTLtotal=length(gr_qtl),
                         QTLs_in_peak = length(unique(df.overlap$snp)), 
                         Pekas_w_qtl = nb_peak_w_qtl, 
                         Peaks_wo_qtl = nb_peak_wo_qtl)
        dt_overlap <- rbind(dt_overlap, dt)
      }
    }
  }
  return(dt_overlap)
}

# Fisher.test
run_fisherTest <- function(df, peakA, peakB){
  fisher_out <- NULL
  # make contingency table -
  #if row#1 vs row#2 has a OR on col#1 (given col1, row1 has a OR when compared with row2)
  label <- paste(peakA, "vs", peakB, sep = " ")
  qtl_cts <- unique(df$QTL_ct)
  for (qtl_ct in qtl_cts){
    dtss <- df[df$QTL_ct %in% qtl_ct, ]
    rownames(dtss) <- dtss$peakT
    dynamicP_w_qtl <- dtss[peakA,"Pekas_w_qtl"]
    dynamicP_wo_qtl <- dtss[peakA,"Peaks_wo_qtl"]
    invarP_w_qtl <- dtss[peakB,"Pekas_w_qtl"]
    invarP_wo_qtl <- dtss[peakB,"Peaks_wo_qtl"]
    contingency_table <- matrix( c(dynamicP_w_qtl, dynamicP_wo_qtl,
                                   invarP_w_qtl,invarP_wo_qtl), nrow=2,
                                 dimnames = list(PeakType = c(peakA, peakB),
                                                 QTL = c("With","Without")))
    print(contingency_table)
    fisher_test <- fisher.test(contingency_table)
    
    dtout <- data.frame(atac_ct=ct_match[qtl_ct,]$atac_ct,
                        qtl_ct=qtl_ct, label,
                        OR=fisher_test$estimate, 
                        pval=fisher_test$p.value)
    fisher_out <- rbind(fisher_out, dtout)
  }
  rownames(fisher_out) <- NULL
  return(fisher_out)
}

qlt_x_peak_sum <- NULL
fisher_sum <- NULL # tidy up all output from diff QTL resources

#---eQTL_Malhotra
dt_used <- eQTL_Malhotra_w_pos[,c("snp","snp_chr","snp_pos_hg38","celltype","eGene")]
colnames(dt_used) <- c("snp","chr","start","celltype","gene")
dt_used <- dt_used[!is.na(dt_used$chr), ]
dt_peak_wQTL <- count_Peak_wQTL(dt_used, "eQTL_Malhotra") 
dt_peak_wQTL[dt_peak_wQTL$Pekas_w_qtl == 0,]
dt_peak_wQTL[dt_peak_wQTL$peakT %in% "dynamicP",]
unique(dt_peak_wQTL$peakT)

#fishe test
fisher_CRE <- run_fisherTest(dt_peak_wQTL, "CRE", "Distal")
fisher_dy <- run_fisherTest(dt_peak_wQTL, "dynamicP", "invariantP")
fisher_dyEnh <- run_fisherTest(dt_peak_wQTL, "Enhancer.dynamicP", "Enhancer.invariantP")
fisher_dy_woDistal <- run_fisherTest(dt_peak_wQTL, "wo_distal.dynamicP", "wo_distal.invariantP")
fisher_dy_dyCRE <- run_fisherTest(dt_peak_wQTL, "CRE.dynamicP", "CRE.invariantP")

fisher_out<- rbind(fisher_CRE, fisher_dy,fisher_dyEnh, fisher_dy_dyCRE,fisher_dy_woDistal)
fisher_out$QTLresource <- rep("sneQTL_Malhotra",dim(fisher_out)[1])

fisher_sum <- rbind(fisher_out, fisher_sum)

#plot
fisher_sum <- read.csv(file.path(outp, "fishertest.QTLs_in_peakset_GeneMatch.csv"),
                       row.names = 1)
dt_fisher <- fisher_sum[fisher_sum$QTLresource %in% "sneQTL_Malhotra",]
unique(dt_fisher$label)
update_label <- data.frame(
  raw=unique(dt_fisher$label)[-1],
  new=c("All Peaks",
        "Enhancer",
        "Enhancer + Promoter",
        "Enhancer + Promoter + Intronic + Exonic")
)
rownames(update_label) <- update_label$raw
dt_fisher$type <- update_label[dt_fisher$label,]$new
dt_fisher <- dt_fisher[!is.na(dt_fisher$type),]

dt_fisher$p_flag <- ifelse(dt_fisher$pval < 0.05, "*", NA)

dt_fisher$type <- factor(dt_fisher$type, 
                         levels = c("All Peaks",
                                    "Enhancer + Promoter + Intronic + Exonic",
                                  "Enhancer + Promoter","Enhancer"))
p <- ggplot(dt_fisher, aes(x = OR, y = -log10(pval))) + 
  geom_point(aes(color = atac_ct, shape=type), size=3) +
  labs(title="QTL Enrichment in Dynamic/Invariant Peaks 
       (Selecting QTLs by matching eGenes to peakGenes)", 
       x="Odds Ratio",y="Fisher Exact p-value (-log10)")+
  xlim(c(0.5,1.5)) +
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
pdf(file.path(outp, "enrich.sneQTL_Malhotra_fisherTest.pdf"), width = 6, height = 4)
p
dev.off()

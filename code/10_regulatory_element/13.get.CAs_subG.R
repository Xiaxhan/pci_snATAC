library(ArchR)
library(parallel)
addArchRThreads(threads = 10)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

add<-theme_set(theme_bw())+
  theme(axis.title.x=element_text(size=10,face = "bold"),
        axis.title.y=element_text(size=10,vjust = 2,hjust = 0.5,face = "bold"),
        axis.text.x=element_text(size=10, face="bold", angle = 45),
        axis.text.y=element_text(size=10,face="bold"),
        legend.position = "top",
        #legend.position = c(1,1),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        panel.grid =element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=1, colour = "black"))

projdir <- "projrmSubset_subPeak" 
proj <- loadArchRProject(projdir)

#Load in co-accessible peaks (for cell types)
proj$subC_Dx  <- paste(proj$subClusters, proj$Clinical.Dx, sep = "-")
subCs <- unique(sort(proj$subClusters))[-62]
cAloops <- list()
for (subc in subCs){
  cAloops[[subc]] <- readRDS(file.path(getwd(), projdir, "CRE", paste("cAloops_",subc, "_rds",  sep = "")))
}

#1. count loop existState of each group
loci.loops <- NULL
gr.lociloops <- NULL
for (subc in names(cAloops)){
  for (subc_condition in names(cAloops[[subc]])){
    print(subc_condition)
    e_cAloops <-  cAloops[[subc]][[subc_condition]]
    addAnote <- rep(subc_condition, length(e_cAloops$CoAccessibility))
    values(e_cAloops$CoAccessibility) <- DataFrame(cellgroup = addAnote, 
                                                   values(e_cAloops$CoAccessibility))
    
    loci <- paste(seqnames(e_cAloops$CoAccessibility),  ranges(e_cAloops$CoAccessibility), sep = "-")
    if(isEmpty(loci.loops)){
      loci.loops <- loci
      gr.lociloops <- e_cAloops$CoAccessibility
    }else{ 
      loci.loops <- c(loci.loops, loci)
      gr.lociloops <- c(gr.lociloops, e_cAloops$CoAccessibility)
    }
  }
}
loci.loops <- unique(sort(loci.loops))
gr.lociloops <- unique(sort(gr.lociloops))

loci.loops<-readRDS(file.path(getwd(), projdir, "CRE", "lociloops_rds"))
gr.lociloops <- readRDS(file.path(getwd(), projdir, "CRE", "grlociloops_rds"))

dt.loopsExist <- data.frame(loci = loci.loops)
for (subc in names(cAloops)){
  for (subc_condition in names(cAloops[[subc]])){
    #subc_condition <- "opc.C9-bvFTD"
    e_cAloops <-  cAloops[[subc]][[subc_condition]]
    
    rawRowname <- colnames(dt.loopsExist)
    loci <- paste(seqnames(e_cAloops$CoAccessibility),  ranges(e_cAloops$CoAccessibility), sep = "-")
    existState <- loci.loops %in% loci
    existState <- as.integer(as.logical(existState)) #!!!!!!
    dt.loopsExist <- cbind(dt.loopsExist, existState)
    colnames(dt.loopsExist) <- c(rawRowname, subc_condition)
  }
}

rownames(dt.loopsExist) <- dt.loopsExist$loci
dt.loopsExist <- dt.loopsExist[,-1]

#2 loops classification: (1) one promoter (2) two promoters (3) none
# fetch both loci for one loop -> loci in a promoter-Peak? -> determine loop type
loopStart <- GRanges(
  seqnames = Rle(seqnames(gr.lociloops)),
  ranges = IRanges(ranges(gr.lociloops)@start)
)
loopEnd <- GRanges(
  seqnames = Rle(seqnames(gr.lociloops)),
  ranges = IRanges(ranges(gr.lociloops)@start + ranges(gr.lociloops)@width - 1)
)  

peaks_W_P <- proj@peakSet[proj@peakSet$peakType %in% "Promoter",]

unique.loopStart <- unique(loopStart)
unique.loopEnd <- unique(loopEnd)

o1 <- findOverlaps(unique.loopStart, peaks_W_P)
loopStart_W_P <- unique.loopStart[queryHits(o1)] 
length(loopStart_W_P) #43514

o2 <- findOverlaps(unique.loopEnd, peaks_W_P)
loopEnd_W_P <- unique.loopEnd[queryHits(o2)] # loopEnd_W_P 43525
length(loopEnd_W_P)


# loop type (In 0/1/2 promoter-peak)
dt.looptype <- NULL
dt.looptype <- data.frame(loci =  paste(seqnames(loopStart), loopStart@ranges, loopEnd@ranges, sep = "-"),
                          leftPoint = paste(seqnames(loopStart), loopStart@ranges, sep = "-"),
                          rightPoint = paste(seqnames(loopEnd), loopEnd@ranges, sep = "-"))
head(dt.looptype)

#Point-peakRange
#L
po1 <- findOverlaps(unique.loopStart, proj@peakSet)
loopStart_P <- unique.loopStart[queryHits(po1)] 

tmpL<- proj@peakSet[subjectHits(po1)]
left.peaks <- data.frame(point = paste(seqnames(loopStart_P), loopStart_P@ranges, sep = "-"), 
                         peak = paste(seqnames(tmpL), tmpL@ranges, sep = "-"),
                         nearestGene = tmpL$nearestGene)

rownames(left.peaks) <- left.peaks$point
left.peaksAll <- left.peaks[dt.looptype$leftPoint,]

#R
po2 <- findOverlaps(unique.loopEnd, proj@peakSet)
loopEnd_P <- unique.loopEnd[queryHits(po2)] 
length(loopEnd_P)

tmpR<- proj@peakSet[subjectHits(po2)]

right.peaks <- data.frame(point = paste(seqnames(loopEnd_P), loopEnd_P@ranges, sep = "-"),
                          peak = paste(seqnames(tmpR), tmpR@ranges, sep = "-"),
                          nearestGene = tmpR$nearestGene)
rownames(right.peaks) <- right.peaks$point
dim(right.peaks)
right.peaksAll <- right.peaks[dt.looptype$rightPoint,]
dim(right.peaksAll)
dt.looptype <- cbind(dt.looptype, data.frame(left_peak=left.peaksAll$peak, 
                                             left_nearestGene = left.peaksAll$nearestGene,
                                             right_peak=right.peaksAll$peak,
                                             right_nearestGene = right.peaksAll$nearestGene))

#Flag W-Promoter
v.startWP <-  paste(seqnames(loopStart_W_P), loopStart_W_P@ranges, sep = "-")
v.endWP <-  paste(seqnames(loopEnd_W_P), loopEnd_W_P@ranges, sep = "-")
dt.looptype <- cbind(dt.looptype, left_W_Promoter = as.integer(as.logical(dt.looptype$leftPoint %in% v.startWP)),
                     right_W_Promoteer = as.integer(as.logical(dt.looptype$rightPoint %in% v.endWP)) )
head(dt.looptype)
#dt.looptype
dt.looptype <- cbind(dt.looptype, nb_ppeak = rowSums(dt.looptype[,8:9]))
head(dt.looptype)

#- get 'with-1pp' as candidate loops
pp1.cAloop <- dt.looptype[dt.looptype$nb_ppeak == 1,]

saveRDS(pp1.cAloop, file.path(getwd(), projdir, "CRE", "subC.cAs_1PromoterPeak.rds"))


#3 plot distriutions of pre candidates enhancers (cAloops with one promoter)
head(dt.loopsExist)

pp1.loopsExist <- dt.loopsExist[rownames(dt.loopsExist) %in% pp1.loci,]

### barplots and heatmaps
pdt <- table(rowSums(pp1.loopsExist))
pdt <- as.data.frame(pdt)
pdt
p <- ggplot(data = pdt, aes(x=Var1, y = Freq)) + geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Number of loops", x = "Number of groups share a loop") + 
  geom_text(mapping = aes(label = Freq), size = 2, vjust = -0.2) + add
pdf(file.path(getwd(), projdir, "Plots", "barplot.distr_ppLoops.pdf"), width = 40)
p
dev.off()

tmpdt <- cbind(pp1.loopsExist, Nb_groups = rowSums(pp1.loopsExist))
tmpdt <- tmpdt[tmpdt$Nb_groups > 1,]
tmpdt <- tmpdt[,-dim(tmpdt)[2]]

head(pp1.loopsExist)
head(as.matrix(pp1.loopsExist))
mt <- as.matrix(pp1.loopsExist)

mt <- as.matrix(tmpdt)

g <- pheatmap(mt, 
              legend=TRUE,
              breaks =  c(0, 1), color =  c("white","red"),
              #legend_breaks = c(0, 1), legend_labels = c(0, 1),
              fontsize = 8,fontsize_row = 8,fontsize_col = 10,
              cluster_cols = T, cluster_rows = T,
              show_rownames = F, show_colnames = T, 
              scale="none", 
              clustering_method='ward.D2', main=" ",
              #display_numbers = TRUE, fontsize_number = 9,
              border_color = NA)

library(ArchR)
library(parallel)
addArchRThreads(threads = 20)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

projdir <- "projrmSubset_subPeak" 
proj <- loadArchRProject(projdir)

#Load in co-accessible peaks (for cell types)
proj$subC_Dx  <- paste(proj$subClusters, proj$Clinical.Dx, sep = "-")
subCs <- unique(sort(proj$subClusters))[-62]
P2Gloops <- list()
for (subc in subCs){
  P2Gloops[[subc]] <- readRDS(file.path(getwd(), projdir, "CRE", paste("P2Gloops_",subc, "_rds",  sep = "")))
}

#0 count loop existState of each group
p2g.loci.loops <- NULL
p2g.gr.lociloops <- NULL
for (subc in names(P2Gloops)){
  for (subc_condition in names(P2Gloops[[subc]])){
    #subc_condition <- "opc.C9-bvFTD"
    print(subc_condition)
    e_p2Gloops <-  P2Gloops[[subc]][[subc_condition]]
    
    loci <- paste(seqnames(e_p2Gloops$Peak2GeneLinks),  ranges(e_p2Gloops$Peak2GeneLinks), sep = "-")
    if(isEmpty(p2g.loci.loops)){
      p2g.loci.loops <- loci
      p2g.gr.lociloops <- e_p2Gloops$Peak2GeneLinks
    }else{ 
      p2g.loci.loops <- c(p2g.loci.loops, loci)
      p2g.gr.lociloops <- c(p2g.gr.lociloops, e_p2Gloops$Peak2GeneLinks)
    }
  }
}
p2g.loci.loops <- unique(sort(p2g.loci.loops))
p2g.gr.lociloops <- unique(sort(p2g.gr.lociloops))

#1 assign loci to gene or peak
peaks_W_P <- proj@peakSet[proj@peakSet$peakType %in% "Promoter",]
gr.genes <- proj@geneAnnotation$genes

loopStart <- GRanges(
  seqnames = Rle(seqnames(p2g.gr.lociloops)),
  ranges = IRanges(ranges(p2g.gr.lociloops)@start, ranges(p2g.gr.lociloops)@start))
loopEnd <- GRanges(
  seqnames = Rle(seqnames(p2g.gr.lociloops)),
  ranges = IRanges(ranges(p2g.gr.lociloops)@start + ranges(p2g.gr.lociloops)@width - 1,
                   ranges(p2g.gr.lociloops)@start + ranges(p2g.gr.lociloops)@width - 1)
)  

#[The peaks] for left-right point in a loop
#left 
op1 <- findOverlaps(loopStart, proj@peakSet)
left_W_peaks <- loopStart[queryHits(op1)]

left_peaksRange <- proj@peakSet[subjectHits(op1)]

left.point2peak <- data.frame(point = paste(seqnames(left_W_peaks), left_W_peaks@ranges, sep = "-"),
                              peak = paste(seqnames(left_peaksRange), left_peaksRange@ranges, sep = "-"),
                              peak_nearestGene = left_peaksRange$nearestGene,
                              peak_type = left_peaksRange$peakType)
head(left.point2peak)
left_WO_peaks <- loopStart[setdiff(1:length(loopStart),queryHits(op1))]
left.point2peak <- rbind(left.point2peak,
                       data.frame(point = paste(seqnames(left_WO_peaks), 
                                                     left_WO_peaks@ranges, sep = "-"),
                                  peak = rep("/", length(left_WO_peaks)),
                                  peak_nearestGene = rep("/", length(left_WO_peaks)),
                                  peak_type = rep("/", length(left_WO_peaks))))
dim(left.point2peak)
head(left.point2peak)

#right
op2 <- findOverlaps(loopEnd, proj@peakSet)

right_W_peaks <- loopEnd[queryHits(op2)]
right_peaksRange <- proj@peakSet[subjectHits(op2)]

right.point2peak <- data.frame(point = paste(seqnames(right_W_peaks), right_W_peaks@ranges, sep = "-"),
                              peak = paste(seqnames(right_peaksRange), right_peaksRange@ranges, sep = "-"),
                              peak_nearestGene = right_peaksRange$nearestGene,
                              peak_type = right_peaksRange$peakType)
head(right.point2peak)
right_WO_peaks <- loopEnd[setdiff(1:length(loopEnd),queryHits(op2))]
right.point2peak <- rbind(right.point2peak,
                         data.frame(point = paste(seqnames(right_WO_peaks), 
                                                       right_WO_peaks@ranges, sep = "-"),
                                    peak = rep("/", length(right_WO_peaks)),
                                    peak_nearestGene = rep("/", length(right_WO_peaks)),
                                    peak_type = rep("/", length(right_WO_peaks))
                         ))
dim(right.point2peak)

#[The genes] for left-right point in a loop
og1 <- findOverlaps(loopStart, gr.genes)

tmpLeft <- unique(loopStart)
otmp <- findOverlaps(tmpLeft, gr.genes)
query.matchMultiGenes <- queryHits(otmp)[duplicated(queryHits(otmp))]
query.matchMultiGenes <- unique(sort(query.matchMultiGenes))
query.matchOneGenes <- setdiff(unique(queryHits(otmp)), query.matchMultiGenes)

# Points WO-G
query.WO_G <- setdiff(1:length(tmpLeft), unique(queryHits(otmp)))

left.point_WO_G <- paste(seqnames(tmpLeft[query.WO_G] ),
                         tmpLeft[query.WO_G] @ranges, sep = "-")

left.point2MultiGene <- NULL
for (qid in query.matchMultiGenes){
  #qid <- 6
  o <- otmp[queryHits(otmp) %in% qid]
  point <- unique(tmpLeft[queryHits(o)])
  point <- paste(seqnames(point), point@ranges, sep = "-")

  matchGenes <- paste(gr.genes[subjectHits(o)]$symbol, collapse  = ";")
  left.point2MultiGene <- rbind(left.point2MultiGene, 
                                data.frame(point = point,
                                           matchGenes = matchGenes))
}
left.point2MultiGene

left.point2uniqueGene <- data.frame(point = paste(seqnames(tmpLeft[query.matchOneGenes]),
                                                  tmpLeft[query.matchOneGenes]@ranges, sep = "-"),
                                    matchGenes = gr.genes[subjectHits(otmp[queryHits(otmp) %in% query.matchOneGenes])]$symbol)

#right
og2 <- findOverlaps(loopEnd, gr.genes)

tmpRight <- unique(loopEnd) 
o.tmpL <- findOverlaps(tmpRight, gr.genes) 
Rquery.matchMultiGenes <- queryHits(o.tmpL)[duplicated(queryHits(o.tmpL))]
Rquery.matchMultiGenes <- unique(sort(Rquery.matchMultiGenes))
Rquery.matchOneGenes <- setdiff(unique(queryHits(o.tmpL)), Rquery.matchMultiGenes)

# Points WO-G
Rquery.WO_G <- setdiff(1:length(tmpRight), unique(queryHits(o.tmpL)))
right.point_WO_G <- paste(seqnames(tmpRight[Rquery.WO_G] ),
                          tmpRight[Rquery.WO_G]@ranges, sep = "-")

right.point2MultiGene <- NULL
j<-0
for (qid in Rquery.matchMultiGenes){ #qid <- 6
  o <- o.tmpL[queryHits(o.tmpL) %in% qid]
  point <- unique(tmpRight[queryHits(o)])
  point <- paste(seqnames(point), point@ranges, sep = "-")
  
  j <- j+1
  print(j)
  print("----")
  print(qid)
  
  matchGenes <- paste(gr.genes[subjectHits(o)]$symbol, collapse  = ";")
  right.point2MultiGene <- rbind(right.point2MultiGene, 
                                data.frame(point = point,
                                           matchGenes = matchGenes))
}
right.point2uniqueGene <- data.frame(point = paste(seqnames(tmpRight[Rquery.matchOneGenes]),
                                                  tmpRight[Rquery.matchOneGenes]@ranges, sep = "-"),
                                    matchGenes = gr.genes[subjectHits(o.tmpL[queryHits(o.tmpL) %in% Rquery.matchOneGenes])]$symbol)

# Combine in a dt
# (1)links; (2)left point (3)left_W_peak?(4) peaktype(5) peakNearestGene (6) left_W_gene?(7) Gene Name? (8)-(13) Right point
dt.P2Ginfo <- data.frame(links = paste(seqnames(loopStart), loopStart@ranges, loopEnd@ranges, sep = "-"),
                         leftPoint = paste(seqnames(loopStart), loopStart@ranges, sep = "-"),
                         rightPoint = paste(seqnames(loopEnd), loopEnd@ranges, sep = "-"))

left.point2peak <- unique(left.point2peak)
rownames(left.point2peak) <- left.point2peak$point
tmpAdd.leftPeak <- left.point2peak[dt.P2Ginfo$leftPoint,]

right.point2peak <- unique(right.point2peak)
rownames(right.point2peak) <- right.point2peak$point
tmpAdd.rightPeak <- right.point2peak[dt.P2Ginfo$rightPoint,]

dt.P2Ginfo <- cbind(dt.P2Ginfo, data.frame(left_W_peak = tmpAdd.leftPeak$peak,
                                           left_peakType = tmpAdd.leftPeak$peak_type,
                                           left_peakNearestGene = tmpAdd.leftPeak$peak_nearestGene,
                                           right_W_peak = tmpAdd.rightPeak$peak,
                                           right_peakType = tmpAdd.rightPeak$peak_type,
                                           right_peakNearestGene = tmpAdd.rightPeak$peak_nearestGene))
left.point2Gsymbol <- NULL
left.point2Gsymbol <-  data.frame(point = c(left.point2uniqueGene$point, left.point2MultiGene$point, left.point_WO_G),
                                  matchGenes= c(left.point2uniqueGene$matchGenes, left.point2MultiGene$matchGenes,
                                                rep("/", length(left.point_WO_G))))

rownames(left.point2Gsymbol) <- left.point2Gsymbol$point

right.point2Gsymbol <- NULL
right.point2Gsymbol <-  data.frame(point = c(right.point2uniqueGene$point, right.point2MultiGene$point, right.point_WO_G),
                                  matchGenes= c(right.point2uniqueGene$matchGenes, right.point2MultiGene$matchGenes,
                                                rep("/", length(right.point_WO_G))))

rownames(right.point2Gsymbol) <- right.point2Gsymbol$point

tmpAdd.leftG <- left.point2Gsymbol[dt.P2Ginfo$leftPoint,]
tmpAdd.rightG <- right.point2Gsymbol[dt.P2Ginfo$rightPoint,]

head(tmpAdd.leftG)
head(tmpAdd.rightG)
dim(tmpAdd.rightG)

dt.P2Ginfo <- cbind(dt.P2Ginfo, data.frame(left_W_gene = tmpAdd.leftG$matchGenes,
                                           right_W_gene = tmpAdd.rightG$matchGenes))

saveRDS(dt.P2Ginfo, file.path(getwd(), projdir, "CRE","subC.P2Ginfo.rds"))


#-------------- Loop Classification: select real "1promoter-1gene"
dt.P2Ginfo <- readRDS(file.path(getwd(), projdir, "CRE","subC.P2Ginfo.rds"))

nb_WO_G <- as.integer(dt.P2Ginfo$left_W_gene %in% "/") + as.integer(dt.P2Ginfo$right_W_gene %in% "/")
dt.P2Ginfo <- cbind(dt.P2Ginfo, L_WO_G = as.integer(dt.P2Ginfo$left_W_gene %in% "/"), 
                    R_WO_G = as.integer(dt.P2Ginfo$right_W_gene %in% "/"), nb_WO_G = nb_WO_G)
dim(dt.P2Ginfo)
head(dt.P2Ginfo)


dim(dt.P2Ginfo[dt.P2Ginfo$left_peakType %in% "Promoter",])

#(1) Althoug a peak is in targetGene, it is still promoter_peak : InGene = nearestGene = linkedGene
dt.peakInGene <- dt.P2Ginfo[dt.P2Ginfo$nb_WO_G == 0,]
dt.peakInGene <- dt.peakInGene[dt.peakInGene$left_peakType %in% "Promoter",]
kept.peak_gene1 <- NULL
for (i in 1:dim(dt.peakInGene)[1]){ #length(p2g.peak_gene)
  if( (dt.peakInGene[i,]$left_peakNearestGene %in% dt.peakInGene[i,]$left_W_gene) & 
      (dt.peakInGene[i,]$left_peakNearestGene %in% dt.peakInGene[i,]$right_W_gene) ){
    
    line <- dt.peakInGene[i,]
    
    if(isEmpty(kept.peak_gene1)){
      kept.peak_gene1 <- line
    }else{
      kept.peak_gene1 <- rbind(kept.peak_gene1,line)}
  }
}

## (2) (A peak is not in gene) promoter-peak's NeaarestGene=linkedGene
dt.peakNotInGene <- dt.P2Ginfo[dt.P2Ginfo$nb_WO_G == 1,]
dt.peakNotInGene <- dt.peakNotInGene[dt.peakNotInGene$L_WO_G == 1 & dt.peakNotInGene$left_peakType %in% "Promoter",]

kept.peak_gene2 <- NULL
for (i in 1:dim(dt.peakNotInGene)[1]){ #length(dt.peakNotInGene)
  if(dt.peakNotInGene[i,]$left_peakNearestGene %in% dt.peakNotInGene[i,]$right_W_gene){
    
    line <- dt.peakNotInGene[i,]
    
    if(isEmpty(kept.peak_gene2)){
      kept.peak_gene2 <- line
    }else{
        kept.peak_gene2 <- rbind(kept.peak_gene2,line)}
  }
}

kept.peak_gene <- rbind(kept.peak_gene1, kept.peak_gene2)
saveRDS(kept.peak_gene, file.path(getwd(), projdir, "CRE","subC.keptP2G.peak2gene.rds"))


######## - [The peaks]
o1 <- findOverlaps(loopStart, peaks_W_P)
loopStart_W_P <- loopStart[queryHits(o1)]
o2 <- findOverlaps(loopEnd, peaks_W_P)
loopEnd_W_P <- loopEnd[queryHits(o2)]

og1 <- findOverlaps(loopStart, gr.genes)
loopStart_W_G <- loopStart[queryHits(og1)]
og2 <- findOverlaps(loopEnd, gr.genes)
loopEnd_W_G <- loopEnd[queryHits(og2)]

findOverlaps(loopEnd[1,], gr.genes)
gr.genes[6:7,]

## Each locus could be 'W-G' / 'W-P' / None
## Left-right: Gene-Peak; Peak-Gene; ...
# loop type (In 0/1/2 promoter-peak)
dt.p2Glooptype <- NULL
dt.p2Glooptype <- data.frame(loci =  paste(seqnames(loopStart), loopStart@ranges, loopEnd@ranges, sep = "-"),
                            left = paste(seqnames(loopStart), loopStart@ranges, sep = "-"),
                            right = paste(seqnames(loopEnd), loopEnd@ranges, sep = "-"))
head(dt.p2Glooptype)
v.startWP <-  paste(seqnames(loopStart_W_P), loopStart_W_P@ranges, sep = "-")
v.endWP <-  paste(seqnames(loopEnd_W_P), loopEnd_W_P@ranges, sep = "-")
v.startWG <-  paste(seqnames(loopStart_W_G), loopStart_W_G@ranges, sep = "-")
v.endWG <-  paste(seqnames(loopEnd_W_G), loopEnd_W_G@ranges, sep = "-")

dt.p2Glooptype <- cbind(dt.p2Glooptype, 
                        left_WP = as.integer(as.logical(dt.p2Glooptype$left %in% v.startWP)),
                        right_WP = as.integer(as.logical(dt.p2Glooptype$right %in% v.endWP)),
                        left_WG = as.integer(as.logical(dt.p2Glooptype$left %in% v.startWG)),
                        right_WG = as.integer(as.logical(dt.p2Glooptype$right %in% v.endWG)),
                        
                        left_peak =
                        left_peak_nearestGene = peaks_W_P$nearestGene,
                        
                        )


saveRDS(dt.p2Glooptype, file.path(getwd(), projdir, "CRE", "dt_p2Glooptype_rds"))

#------
for (subc in names(P2Gloops)){
  for (subc_condition in names(P2Gloops[[subc]])){
    #subc_condition <- "opc.C9-bvFTD"
    print(subc_condition)
    e_p2Gloops <-  P2Gloops[[subc]][[subc_condition]]
  }
}



#save loop_existence_state for each group: 'dt.loopsExist' (row-loopLoci col-celltypeCondition)
p2g.dt.loopsExist <- data.frame(loci = p2g.loci.loops)
for (subc in names(P2Gloops)){
  for (subc_condition in names(P2Gloops[[subc]])){
    e_P2Gloops <-  P2Gloops[[subc]][[subc_condition]]
    
    rawRowname <- colnames(p2g.dt.loopsExist)
    loci <- paste(seqnames(e_P2Gloops$Peak2GeneLinks),  ranges(e_P2Gloops$Peak2GeneLinks), sep = "-")
    existState <- p2g.loci.loops %in% loci
    existState <- as.integer(as.logical(existState)) #!!!!!!
    p2g.dt.loopsExist <- cbind(p2g.dt.loopsExist, existState)
    colnames(p2g.dt.loopsExist) <- c(rawRowname, subc_condition)
  }
}

rownames(p2g.dt.loopsExist) <- p2g.dt.loopsExist$loci
p2g.dt.loopsExist <- p2g.dt.loopsExist[,-1]
saveRDS(p2g.dt.loopsExist, file.path(getwd(), projdir, "CRE", "P2G_dtloopsExist_rds"))




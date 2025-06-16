# Extract bed files of cluster specific peaks for LDSC
library(GenomicRanges)
library(ArchR)
library(parallel)
addArchRThreads(threads = 20)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")
library(ggplot2)
library(ggpubr)

projdir <- "projrmSubset_subPeak"
proj <- loadArchRProject(projdir)

clusters <- unique(proj$subClusters)

peakSet.cluster <- list()
wholePeakSet <- NULL
nb_peaks <- NULL
for (clu in clusters){
  fname <- paste(clu, "-reproduciblePeaks.gr.rds", sep = "")
  Raw_peaks <- readRDS(file.path(projdir, "PeakCalls", fname))
  m <- findMatches(Raw_peaks, proj@peakSet)
  o <- findOverlaps(Raw_peaks, proj@peakSet)
  MapPeaks <- unique(proj@peakSet[subjectHits(o)])

  nb_matchPeaks <- length(unique(queryHits(m)))
  nb_overlapPeaks <- length(unique(queryHits(o)))
  nb_clusterPeaks <- length(Raw_peaks)
  nb_mapPeaks <- length(MapPeaks)
  nb_peaks <- rbind(nb_peaks, data.frame(clu = clu,
                                         nb_clusterPeaks=nb_clusterPeaks,
                                         nb_matchPeaks=nb_matchPeaks,
                                         nb_overlapPeaks = nb_overlapPeaks,
                                         nb_mapPeaks = nb_mapPeaks
                                         ))
  peaks <- Raw_peaks 
  peakSet.cluster[[clu]] <- Raw_peaks

  values(peaks) <- data.frame(cluster = rep(clu, length(peaks)),
                              pLoca =  paste(seqnames(peaks), ranges(peaks), sep = "-") )
  if(isEmpty(wholePeakSet)){ 
    wholePeakSet <- peaks
  }else{ 
    wholePeakSet <- c(wholePeakSet, peaks)
    }
}
nb_peaks
sum(nb_peaks$nb_clusterPeaks)

wholePeakSet <- unique(wholePeakSet) 
wholePeakSet


wholePeakSet <- wholePeakSet[ wholePeakSet$cluster %ni% "undefined"]
peakSet.cluster$undefined <- NULL

peakSet.clusterSpecific <- list()
for (clu in sort(names(peakSet.cluster))){
  peaks <- peakSet.cluster[[clu]]
  otherPeakSet <-  wholePeakSet[ wholePeakSet$cluster %ni% clu]
  o <- findOverlaps(peaks, otherPeakSet)
  index.unique <- setdiff(1:length(peaks), unique(queryHits(o)) )
  uniquePeaks <- peaks[index.unique]
  peakSet.clusterSpecific[[clu]] <-uniquePeaks
}

counts <- NULL
for (cluster in sort(names(peakSet.cluster)) ){
  nb_all <- length(peakSet.cluster[[cluster]])
  if(is.null(nb_all)){nb_all <- 0}
  nb <- length(peakSet.clusterSpecific[[cluster]])
  if(is.null(nb)){nb <- 0}
  counts <- rbind(counts, data.frame(cluster=cluster, all=nb_all, unique=nb))
}
counts$cluster <- factor(counts$cluster , levels=sort(counts$cluster) )

outp <- "cluster_specific_atacPeak"
for (clu in names(peakSet.clusterSpecific)){
  keep_peaks <- peakSet.clusterSpecific[clu]
  cluster_specific_peak <- data.frame(keep_peaks)[,1:3]
  name <- paste(clu, "bed", sep = ".")
  write.table(cluster_specific_peak, quote = F, sep = "\t", row.names = F, col.names = F, file.path(outp, name))
}
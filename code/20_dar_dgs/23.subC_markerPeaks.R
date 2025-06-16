library(ArchR)
library(parallel)
addArchRThreads(threads = 12)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

projdir <- "projrmSubset_subPeak" 
proj <- loadArchRProject(projdir)

dir <- "mydir"

markersPeaks <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "PeakMatrix", 
  groupBy = "subClusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersPeaks, file = file.path(dir, "markersPeaks_subC.rds"))

markersGS <- getMarkerFeatures(
  ArchRProj = proj, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "subClusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersGS, file = file.path(dir, "markersGeneScore_subC.rds"))

# annotation
markersPeaks <- readRDS(file = file.path(dir, "markersPeaks_subC.rds"))
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.1", returnGR = TRUE)

for (group in names(markerList)){
  overlaps <- findOverlaps(markerList[[group]], proj@peakSet)
  values(markerList[[group]]) <- DataFrame(values(markerList[[group]]),
                                           nearestGene = proj@peakSet[subjectHits(overlaps)]$nearestGene,
                                           peakType = proj@peakSet[subjectHits(overlaps)]$peakType)
  fname <- paste(group, "markerPeaks_fdr0.1.csv", sep = "")
  write.csv(data.frame(markerList[[group]]), file.path(getwd(), projdir, "subC_markers", fname))
}

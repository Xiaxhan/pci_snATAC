# Call marker peaks for disaeses within each subC

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
Dxs <- c("Control","AD","bvFTD","PSP_S")
majorCs <- c("ast","mg","EX","IN","odc","opc")
peakTypes <- c("Distal","Enhancer","Exonic","Intronic","Promoter")
TEclasses <- c("DNA","LINE","LTR","SINE","otherTEs")

colors.ct <- c(ArchRPalettes$stallion[c(1,5,3,11)], "#b9a802","#308ac4")
names(colors.ct) <-c("ast","mg","IN","EX","odc","opc")

projdir <- "projrmSubset_subPeak" 
proj <- loadArchRProject(projdir)
dropSubCs <- c("undefined","ast.C6","mg.C8","mg.C10","mg.C15",
               "mg.C1", "mg.C5", "neu.C3", "neu.C4", "odc.C2", "neu.C2", "mg.C2")
proj$subClusters[proj$subClusters %like% "neuron"] <- gsub("neuron","neu",
                                                           proj$subClusters[proj$subClusters %like% "neuron"])
projclean <- proj[!proj$subClusters %in% dropSubCs,]

projclean$new_majorC <- projclean$subClusters
projclean$new_majorC[projclean$new_majorC %in% 
                       c("neu.C6","neu.C7","neu.C8","neu.C9")] <- "IN" 
projclean$new_majorC[projclean$new_majorC %like% "neu"] <- "EX" 
projclean$new_majorC <- gsub("\\..*","", projclean$new_majorC)
projclean <- projclean[!projclean$Sample %in% c("I1_7","P1_7_at1_7"),]

MarkerPeaks_Raw <- list()
MarkerPeaks_byPval <- list()
for (cluster in clusters){
  print(cluster)
  projss <- projclean[projclean$subClusters %in% cluster,]
  markersPeaks <- getMarkerFeatures(
    ArchRProj = projss, 
    useMatrix = "PeakMatrix", 
    groupBy = "Clinical.Dx",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  MarkerPeaks_Raw[[cluster]] <- markersPeaks
  
  markerList <- getMarkers(markersPeaks, cutOff = "Pval < 0.05", returnGR = TRUE)
  MarkerPeaks_byPval[[cluster]] <- markerList
}
library(ArchR)
library(parallel)
addArchRThreads(threads = 1)
addArchRGenome("hg38")

# load the main object
projss1 <- loadArchRProject("projrmSubset_subPeak")

celltype <- "MG" # or 'AST' 'ODC' 'OPC'
  
cellKeep <- rownames(projss1@cellColData[projss1@cellColData$Clusters_r1_iLSI6 %in% celltype,])
length(cellKeep) / length(rownames(projss1@cellColData)) #0.4078475

# subset and save
subsetArchRProject(
  ArchRProj = projss1,
  cells = cellKeep,
  outputDirectory = "proj_microglia",
  dropCells = TRUE,
  logFile = NULL,
  threads = getArchRThreads(),
  force = TRUE
)

# load subset of a single cell type
proj <- loadArchRProject("proj_microglia")

# run LSI, harmony and UMAP
# modify specifically for each cell type, show example here
lsi_reDim <- "iLSI6" 
iters <- 6
har_reDim <- "cs_harmony_iLSI6"
cluster_reDim <- "Clusters_r2_cs_harmony_iLSI6"
clu_res <- 0.2

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", dimsToUse = 1:30,
                        name = lsi_reDim,
                        varFeatures = 15000,
                        iterations = iters,
                        clusterParams = list(resolution = c(0.1, 0.2, 0.3, 0.4, 0.6),
                                             sampleCells = 10000, n.start = 10)
)

proj <- addHarmony(ArchRProj= proj,
                   reducedDims = lsi_reDim,
                   name = har_reDim,
                   groupBy="PrepBatch") 

proj <- addClusters(input = proj,
                    reducedDims = har_reDim,
                    name = cluster_reDim,
                    resolution = clu_res,
                    force = TRUE)

saveArchRProject(ArchRProj = proj, outputDirectory = projdir, load=F)


# Markers
markersGS <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "GeneScoreMatrix",
  groupBy = cluster_reDim,
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(markersGS, file.path(getwd(), projdir, "marker_genescores.rds"))

#select top 5 (log2FC) unique markers to label
markersGS <- readRDS(file.path(getwd(), projdir, "marker_genescores.rds"))
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

markerGenes <- NULL
for (cid in 1:length(unique(proj$Clusters_r2_cs_harmony_iLSI6))){
  cid <- paste("C", cid, sep = "")
  print(cid)
  #markerList$C1[order(-markerList$C1$Log2FC),][1:5,]
  cid_markers <- markerList[[cid]][order(-markerList[[cid]]$Log2FC),]$name
  markerGenes <- c(markerGenes, cid_markers[1:5])
}
length(markerGenes)
markerGenes <- markerGenes[!is.na(markerGenes)]

heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE,
  clusterCols = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores_marker_Heatmap", width = 6, height = 8, ArchRProj = proj, addDOC = FALSE)
        

# enrichment for marker
library(enrichR)
dbs <- listEnrichrDbs()
usedbs <- dbs$libraryName
dt.enriched <- list()
dt.enriched.clean <- list()
for (ct in c("ast","mg","odc","opc","neu")){
  markersGS <- readRDS(file.path(getwd(), projdir, ct, "marker_genescores.rds"))
  mk <- getMarkers(markersGS, cutOff = "FDR <= 0.1 & Log2FC > 1", returnGR = F)
  for (subC in names(mk) ){
    print(subC)
    #all cells
    dt <- mk[[subC]]
    if(dim(dt)[1]!=0){
      enriches <- enrichr(as.character(dt$name), usedbs)
      #outs
      for (db in usedbs){
        print(db)
        enriched.db <- enriches[[db]]
        enriched.db <- enriched.db[enriched.db$Adjusted.P.value < 0.1 &
                                     enriched.db$P.value < 0.05, ]
        dt.enriched[[ct]] <- rbind(dt.enriched[[ct]],
                                   data.frame(subC = rep(subC, dim(enriched.db)[1]),
                                              db = rep(db, dim(enriched.db)[1]),
                                              enriched.db))
      }
    }
  }
  write.csv(dt.enriched[[ct]], file.path(outp,paste(ct,"enrichR.markers_SubC.txt",sep = ".")))
}
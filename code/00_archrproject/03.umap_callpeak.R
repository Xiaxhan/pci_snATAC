library(ArchR)
library(parallel)
addArchRThreads(threads = 20)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

projdir <- "projALL_rmDoublets"
proj <- loadArchRProject(projdir)

# add LSI
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", dimsToUse = 1:30,
                        name = "iLSI6",
                        varFeatures = 25000,
                        iterations = 6,
                        clusterParams = list(resolution = c(0.1, 0.2, 0.4, 0.6, 0.8),
                                             sampleCells = 10000, n.start = 10)
)

# Hamonry
proj <- addHarmony(ArchRProj= proj,
                   reducedDims = "iLSI6",
                   name = "harmony_iLSI6",
                   groupBy="PrepBatch")

# add UMAP
proj <- addUMAP(nNeighbors = 30,
                minDist = 0.1,
                metric = "cosine",
                ArchRProj = proj,
                reducedDims = "harmony_iLSI6",
                name = "UMAP2_iLSI6",
                force = TRUE
)

# clustering
proj <- addClusters(input = proj, reducedDims = "iLSI6", name = "Clusters_r1_iLSI6",
                    resolution = 0.1, force = TRUE)

#make pseudo-bulk replicates for clusters ( ~min)
proj <- addGroupCoverages(ArchRProj = proj, groupBy = "subClusters")

pathToMacs2 <- "/home/xhan/anaconda3/envs/macs2/bin/macs2"
proj <- addReproduciblePeakSet(
  ArchRProj = proj, groupBy = "subClusters", pathToMacs2 = pathToMacs2)

projnew <- addPeakMatrix(proj)

saveArchRProject(ArchRProj = projnew, outputDirectory = "projrmSubset_subPeak", load=F)
# Call DARs for individual cell types based on the downsampled dataset.

library(ArchR)
library(parallel)
addArchRThreads(threads = 5)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

proj <- loadArchRProject(projdir)

outp <- "mydir"

dropSubCs <- c("undefined","ast.C6","mg.C8","mg.C10","mg.C15",
               "mg.C1", "mg.C5", "neu.C3", "neu.C4", "odc.C2", "neu.C2", "mg.C2")
proj$subClusters[proj$subClusters %like% "neuron"] <- gsub("neuron","neu",
                                                           proj$subClusters[proj$subClusters %like% "neuron"])
projclean <- proj[!proj$subClusters %in% dropSubCs,]
length(unique(projclean$subClusters))

subClusters <- sort(unique(projclean$subClusters))

tmp.majorC <- projclean$subClusters
tmp.majorC[tmp.majorC %in% c("neu.C6","neu.C7","neu.C8", "neu.C9")] <- "IN"
tmp.majorC[tmp.majorC %like% "neu"] <- "EX"

tmp.majorC[tmp.majorC %like% "ast"] <- "ASC"
tmp.majorC[tmp.majorC %like% "mg"] <- "MG"
tmp.majorC[tmp.majorC %like% "odc"] <- "ODC"
tmp.majorC[tmp.majorC %like% "opc"] <- "OPC"
projclean@cellColData <- cbind(projclean@cellColData, majorC = tmp.majorC)
marjoCs <- c("ASC", "MG", "EX", "IN", "ODC", "OPC")
proj <- projclean

region <- "midInsula" # midInsula PreCG
proj.re <- proj[proj$region %in% region]
proj.re

majorC_diagnose <- paste(proj.re$Clinical.Dx, proj.re$majorC, sep = "|")
proj.re$majorC_diagnose <- majorC_diagnose

samples <- unique(proj.re$Sample)

#Record how many samples have highQ cells > 15 per cell type
proj.highQ <- proj.re[proj.re@cellColData$TSSEnrichment > 4,]
proj.highQ$dxSample <- paste(proj.highQ$Clinical.Dx, proj.highQ$Sample, sep = "|")
df <- table(proj.highQ$majorC, proj.highQ$dxSample)
write.csv(t(df), file.path(outp, paste(region, "highQ_majorCcell.csv", sep = ".")) )

#Downsample: 30 nuclei per sample for each cell type
cellList <- list()
for (i in 1:10){
  sampled_cells <- NULL
  for (ct in marjoCs){
    sample_cells <- NULL
    for (sp in samples){
      print(paste("---",i, ct, sp))
      cells_set <- rownames(proj.re[proj.re@cellColData$Sample %in% sp & 
                             proj.re@cellColData$majorC %in% ct,]@cellColData)
      set.seed(i)
      sampled <- sample(cells_set, 30) 
      sampled_cells <- c(sampled_cells, sampled)
    }
  }
  print(paste(i, length(sampled_cells)))
  cellList[[paste("seed",i,sep = "")]] <- sampled_cells
}

#call DAR
diagnoses <- c("AD", "bvFTD", "PSP_S")
majorCell <- c("ASC", "MG", "EX","IN","ODC", "OPC")
for (seed in names(cellList)){
  sampled_cells <- cellList[[seed]]
  proj.sampled <- proj.re[rownames(proj.re@cellColData) %in% sampled_cells]
  print(seed)
  DARtest <- list()
  for (disease in diagnoses){
    for (majortype in majorCell){
      group.disease <- paste(disease, majortype, sep = "|")
      group.control <- paste("Control", majortype, sep = "|")
      print(group.disease)
      
      darTest <- getMarkerFeatures(
        ArchRProj = proj.sampled, 
        useMatrix = "PeakMatrix",
        groupBy = "majorC_diagnose",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = group.disease,
        bgdGroups = group.control
      )
      DARtest[[group.disease]] <- darTest
    }
  }
  
  fileout <-  paste(region, seed, "DAR_6majorC.rds", sep = ".")
  saveRDS(DARtest, file = file.path(outp,  file = fileout ) )
  print("Done!")
}



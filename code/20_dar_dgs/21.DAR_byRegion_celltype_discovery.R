#Split the data into 6 cell types (EX, IN subtypes); 
#perform differential accessibility (DA) analysis for each cell type, stratified by brain region.

library(ArchR)
library(parallel)
addArchRThreads(threads = 10)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

projdir <- "projrmSubset_subPeak" 
proj <- loadArchRProject(projdir)

# clean subtypes 
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

proj <- projclean

# pick region
region <- "midInsula" # midInsula PreCG
proj.re <- proj[proj$region %in% region]
proj.re

majorC_diagnose <- paste(proj.re$Clinical.Dx, proj.re$majorC, sep = "|")
proj.re$majorC_diagnose <- majorC_diagnose


majorCell <- c("ASC", "MG", "ODC", "OPC", "EX", "IN")

#call DAR
for (disease in diagnoses){
  for (majortype in majorCell){

    group.disease <- paste(disease, majortype, sep = "|")
    group.control <- paste("Control", majortype, sep = "|")
    print(group.disease)
    
    darTest <- getMarkerFeatures(
      ArchRProj = proj.re, 
      useMatrix = "PeakMatrix",
      groupBy = "majorC_diagnose",
      testMethod = "wilcoxon",
      bias = c("TSSEnrichment", "log10(nFrags)"),
      useGroups = group.disease,
      bgdGroups = group.control
    )
    groupDx <- paste(disease, majortype, sep = ".")
    fileout <-  paste("DAR_6majorC", region, groupDx, "with50subC.rds", sep = ".")
    saveRDS(darTest, file = file.path(getwd(), projdir,  file = fileout ) )
    print("Done!")
  }
}

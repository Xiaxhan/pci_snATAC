library(ArchR)
library(parallel)
addArchRThreads(threads = 10)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

setwd("/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318")
projdir <- "projrmSubset_subPeak" 
proj <- loadArchRProject(projdir)
 
projmeta <- readRDS("/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318/projATAC_50subC_meta_584904.rds")
finalSubCs <- sort(unique(projmeta$subClusters))
  
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
rm(proj)
majorCs <- sort(unique(projclean$new_majorC))
majorCs
projclean$DX <- projclean$Clinical.Dx
projclean$DX[ projclean$DX %in% "Control" ] <- "1Control"
  
projclean$subC_Dx  <- paste(projclean$subClusters, projclean$Clinical.Dx, sep = "-")

subCs <- unique(sort(projclean$subClusters))[-62]

ctype <- "mg" 
subCs <- finalSubCs[finalSubCs %like% ctype]
subCs

for (subc in subCs){
  print(subc)
  cAs <- list()
  
  ssCells <- projclean$cellNames[projclean$subClusters %in% subc]
  projclean <-  addCoAccessibility(
    ArchRProj = projclean,
    reducedDims = "iLSI6",
    cellsToUse = ssCells
  )
  
    cAs[[subc]] <- getCoAccessibility(
      ArchRProj = projclean,
      corCutOff = 0.5,
      resolution = 1,
      returnLoops = FALSE
    )
  
  for (condi in unique(sort(projclean$Clinical.Dx))){
    subc_condi <- paste(subc, condi, sep = "-")
    print(subc_condi)
    ssCells <- projclean$cellNames[projclean$subC_Dx %in% subc_condi]
    
    if(length(ssCells) < 100){next}
    
    projclean <-  addCoAccessibility(
      ArchRProj = projclean,
      reducedDims = "iLSI6",
      cellsToUse = ssCells
    )
    
      cAs[[subc_condi]] <- getCoAccessibility(
        ArchRProj = projclean,
        corCutOff = 0.5,
        resolution = 1,
        returnLoops = FALSE
      )
    
    
  }
  saveRDS(cAs, file = file.path(getwd(), projdir, "CRE", paste("cAloops_",subc, "_rds",  sep = "")))
  print("Done")
}

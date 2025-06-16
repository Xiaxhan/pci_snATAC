library(ArchR)
library(parallel)
addArchRThreads(threads = 20)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

setwd("/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318")
projdir <- "projrmSubset_subPeak" 
proj <- loadArchRProject(projdir)

if(T){
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
  
}

proj <- projclean

proj$subC_Dx  <- paste(proj$subClusters, proj$Clinical.Dx, sep = "-")


ctype <- "mg"  # !!!!!!!!!!

subCs <- finalSubCs[finalSubCs %like% ctype]

for (subc in subCs){
  print(subc)
  
  P2Gloops <- list()
  
  ssCells <- proj$cellNames[proj$subClusters %in% subc]
  
  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "iLSI6",
    cellsToUse = ssCells
  )
  
  P2Gloops[[subc]] <- getPeak2GeneLinks(
    ArchRProj = proj,
    corCutOff = 0.45,
    resolution = 1,
    FDRCutOff = 0.1,
    returnLoops = TRUE
  )
  
  for (condi in unique(sort(proj$Clinical.Dx))){
    subc_condi <- paste(subc, condi, sep = "-")
    print(subc_condi)
    ssCells <- proj$cellNames[proj$subC_Dx %in% subc_condi]
    
    if(length(ssCells) < 100){next}
    
    proj <- addPeak2GeneLinks(
      ArchRProj = proj,
      reducedDims = "iLSI6",
      cellsToUse = ssCells
    )
    
    P2Gloops[[subc_condi]] <- getPeak2GeneLinks(
      ArchRProj = proj,
      corCutOff = 0.45,
      resolution = 1,
      FDRCutOff = 0.1,
      returnLoops = TRUE
    )
    
  saveRDS(P2Gloops, 
          file = file.path(getwd(), projdir, "CRE/cAs_p2g_cor0.1", paste("P2Gloops_",subc, "_rds",  sep = "")))
  
}


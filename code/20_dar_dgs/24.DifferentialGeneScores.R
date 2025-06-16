library(rhdf5)
library(data.table)
library(dplyr)
library(ArchR)
library(parallel)
library(preprocessCore)
addArchRThreads(threads = 20)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
projdir <- "projrmSubset_subPeak"

projmeta <- readRDS("projATAC_50subC_meta_584904.rds")
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


#subC :
#"DEGtest_atacScore_subCs.rds" 
#"df_dxDEG_fdr0.1_subCs.rds
#"df_dxDEG_Pval0.1_subCs.rds" 
if(T){
  #subCs
  DEGtest <- list()
  for (ct in c("ast","mg", "EX","IN","odc","opc")){
    subCs <- unique(projclean[projclean$new_majorC %like% ct,]$subClusters) %>% sort()
    for (subC in subCs){
      projsubC <- projclean[projclean$subClusters %in% subC,]
      for (dx in c("AD","bvFTD","PSP_S")){
        dmrTest <- getMarkerFeatures(
          ArchRProj = projsubC,
          useMatrix = "GeneScoreMatrix",
          groupBy = "Clinical.Dx",
          testMethod = "wilcoxon",
          bias = c("TSSEnrichment", "log10(nFrags)"),
          useGroups = dx,
          bgdGroups = "Control"
        )
        groupName <- paste(subC, dx, sep = "|")
        DEGtest[[groupName]] <- dmrTest
        print(groupName)
      }
    }
  }
 # saveRDS(DEGtest, file.path(projdir, paste("DEGtest_atacScore_subCs.rds")))
  print("saved")
  
  DEGtest <- readRDS(file.path(projdir, paste("DEGtest_atacScore_subCs.rds")) )
  
  df_differentialGene <- NULL # based on gene scores
  for (groupName in names(DEGtest)){
    dag <- DEGtest[[groupName]]
  #  sigDAG <- getMarkers(dag, cutOff = "FDR <= 0.1", returnGR = F)
    sigDAG <- getMarkers(dag, cutOff = "Pval <= 0.1", returnGR = F)
    dx <- names(sigDAG)
    tmpdt <- sigDAG[[dx]]
    dt_group <- tmpdt[,c("name","Log2FC","FDR","MeanDiff")]
    dt_group <- cbind(dt_group, groupName=rep(groupName, dim(dt_group)[1]))
    df_differentialGene <- rbind(df_differentialGene, dt_group )
  }
  dim(df_differentialGene)
  df_differentialGene$subC <- gsub("\\|.*","",df_differentialGene$groupName)
  df_differentialGene$dx <- gsub(".*\\|","",df_differentialGene$groupName)
#  saveRDS(df_differentialGene, file.path(projdir, paste("df_dxDEG_fdr0.1_subCs.rds")))
  saveRDS(df_differentialGene, file.path(projdir, paste("df_dxDEG_Pval0.1_subCs.rds")))
  head(df_differentialGene)
  
}

#majorC:
#"DEGtest_atacScore_majorCs.rds" 
#"df_dxDEG_fdr0.1_majorCs.rds"
#"df_dxDEG_Pval0.1_majorCs.rds"
if(T){
  DEGtest <- list()
  for (ct in majorCs){
    projct <- projclean[projclean$new_majorC %in% ct,]
    for (dx in c("AD","bvFTD","PSP_S")){
      dmrTest <- getMarkerFeatures(
        ArchRProj = projct,
        useMatrix = "GeneScoreMatrix",
        groupBy = "Clinical.Dx",
        testMethod = "wilcoxon",
        bias = c("TSSEnrichment", "log10(nFrags)"),
        useGroups = dx,
        bgdGroups = "Control"
      )
      groupName <- paste(ct, dx, sep = "|")
      DEGtest[[groupName]] <- dmrTest
      print(groupName)
    }
  }
  #saveRDS(DEGtest, file.path(projdir, paste("DEGtest_atacScore_majorCs.rds")))
  
  DEGtest <- readRDS(file.path(projdir, paste("DEGtest_atacScore_majorCs.rds")) )
  
  df_differentialGene <- NULL # based on gene scores
  for (groupName in names(DEGtest)){
    dag <- DEGtest[[groupName]]
    sigDAG <- getMarkers(dag, cutOff = "FDR <= 0.1", returnGR = F)
  #  sigDAG <- getMarkers(dag, cutOff = "Pval <= 0.1", returnGR = F)
    dx <- names(sigDAG)
    tmpdt <- sigDAG[[dx]]
   # dt_group <- tmpdt[,c("name","Log2FC","FDR")]
    dt_group <- tmpdt[,c("name","Log2FC","FDR","MeanDiff")]
    dt_group <- cbind(dt_group, groupName=rep(groupName, dim(dt_group)[1]))
    df_differentialGene <- rbind(df_differentialGene, dt_group )
  }
  dim(df_differentialGene)
  df_differentialGene$subC <- gsub("\\|.*","",df_differentialGene$groupName)
  df_differentialGene$dx <- gsub(".*\\|","",df_differentialGene$groupName)
  saveRDS(df_differentialGene, file.path(projdir, paste("df_dxDEG_fdr0.1_majorCs.rds")))
#  saveRDS(df_differentialGene, file.path(projdir, paste("df_dxDEG_Pval0.1_majorCs.rds")))
  head(df_differentialGene)
  
}


# by regions - majorC:
#"DEGtest_atacScore_majorCs_region.rds" 
#"df_dxDEG_fdr0.1_majorCs_region.rds"
#"df_dxDEG_Pval0.1_majorCs_region.rds"
if(T){
  DEGtest <- list()
  for (region in c("PreCG","midInsula")){
    projre <- projclean[projclean$region %in% region,]
    for (ct in majorCs){
      projct <- projre[projre$new_majorC %in% ct,]
      for (dx in c("AD","bvFTD","PSP_S")){
        dmrTest <- getMarkerFeatures(
          ArchRProj = projct,
          useMatrix = "GeneScoreMatrix",
          groupBy = "Clinical.Dx",
          testMethod = "wilcoxon",
          bias = c("TSSEnrichment", "log10(nFrags)"),
          useGroups = dx,
          bgdGroups = "Control"
        )
        groupName <- paste(region, ct, dx, sep = "|")
        DEGtest[[groupName]] <- dmrTest
        print(groupName)
      }
    }
  }
  saveRDS(DEGtest, file.path(projdir, paste("DEGtest_atacScore_majorCs_region.rds")))
  
  DEGtest <- readRDS(file.path(projdir, paste("DEGtest_atacScore_majorCs.rds")) )
  
  df_differentialGene <- NULL # based on gene scores
  for (groupName in names(DEGtest)){
    dag <- DEGtest[[groupName]]
   # sigDAG <- getMarkers(dag, cutOff = "FDR <= 0.1", returnGR = F)
    sigDAG <- getMarkers(dag, cutOff = "Pval <= 0.1", returnGR = F)
    dx <- names(sigDAG)
    tmpdt <- sigDAG[[dx]]
    dt_group <- tmpdt[,c("name","Log2FC","FDR","MeanDiff")]
    dt_group <- cbind(dt_group, groupName=rep(groupName, dim(dt_group)[1]))
    df_differentialGene <- rbind(df_differentialGene, dt_group )
  }
  dim(df_differentialGene)
  df_differentialGene$region <- gsub("\\|.*","",df_differentialGene$groupName)
  df_differentialGene$dx <- gsub(".*\\|","",df_differentialGene$groupName)
  df_differentialGene$ct <- gsub(".*\\|(.*?)\\|.*", "\\1",df_differentialGene$groupName)
  #saveRDS(df_differentialGene, file.path(projdir, paste("df_dxDEG_fdr0.1_majorCs_region.rds")))
  saveRDS(df_differentialGene, file.path(projdir, paste("df_dxDEG_Pval0.1_majorCs_region.rds")))
  head(df_differentialGene)
  
}

# by regions - subCs:
#"DEGtest_atacScore_subCs_region.rds" 
#"df_dxDEG_fdr0.1_subCs_region.rds"
# df_dxDEG_Pval0.1_subCs_region.rds
if(T){
  #subCs
  DEGtest <- list()
  for (region in c("PreCG","midInsula")){
    projre <- projclean[projclean$region %in% region,]
    for (ct in majorCs){
      projct <- projre[projre$new_majorC %in% ct,]
      subCs <- unique(projct[projct$new_majorC %like% ct,]$subClusters) %>% sort()
      for (subC in subCs){
        projsubC <- projct[projct$subClusters %in% subC,]
        for (dx in c("AD","bvFTD","PSP_S")){
          dmrTest <- getMarkerFeatures(
            ArchRProj = projsubC,
            useMatrix = "GeneScoreMatrix",
            groupBy = "Clinical.Dx",
            testMethod = "wilcoxon",
            bias = c("TSSEnrichment", "log10(nFrags)"),
            useGroups = dx,
            bgdGroups = "Control"
          )
          groupName <- paste(region, subC, dx, sep = "|")
          DEGtest[[groupName]] <- dmrTest
          print(groupName)
        }
      }
    }
  }
 
  saveRDS(DEGtest, file.path(projdir, paste("DEGtest_atacScore_subCs_region.rds")))
  print("saved")
  
  DEGtest <- readRDS(file.path(projdir, paste("DEGtest_atacScore_subCs_region.rds")) )
  
  df_differentialGene <- NULL # based on gene scores
  for (groupName in names(DEGtest)){
    dag <- DEGtest[[groupName]]
#    sigDAG <- getMarkers(dag, cutOff = "FDR <= 0.1", returnGR = F)
    sigDAG <- getMarkers(dag, cutOff = "Pval <= 0.1", returnGR = F)
    dx <- names(sigDAG)
    tmpdt <- sigDAG[[dx]]
    dt_group <- tmpdt[,c("name","Log2FC","FDR","MeanDiff")]
    dt_group <- cbind(dt_group, groupName=rep(groupName, dim(dt_group)[1]))
    df_differentialGene <- rbind(df_differentialGene, dt_group )
  }
  dim(df_differentialGene)
  df_differentialGene$subC <- gsub("\\|.*","",df_differentialGene$groupName)
  df_differentialGene$dx <- gsub(".*\\|","",df_differentialGene$groupName)
#  saveRDS(df_differentialGene, file.path(projdir, paste("df_dxDEG_fdr0.1_subCs_region.rds")))
  saveRDS(df_differentialGene, file.path(projdir, paste("df_dxDEG_Pval0.1_subCs_region.rds")))
  head(df_differentialGene)
  
}

#output csv

#Differential gene scores
majorC_differentialGene <- readRDS( file.path(projdir, paste("df_dxDEG_fdr0.1_majorCs.rds") ) )
majorC_dGS_region <- readRDS( file.path(projdir, paste("df_dxDEG_fdr0.1_majorCs_region.rds") ) )
cols <- c("name","Log2FC","FDR","groupName","region","ct","dx")
colnames(majorC_differentialGene)[colnames(majorC_differentialGene) %in% "subC"] <- "ct"
majorC_differentialGene$region <- rep("none",dim(majorC_differentialGene)[1])

majorC_DGS <- rbind(majorC_differentialGene[,cols], majorC_dGS_region[,cols])


tmpDT <- as.data.frame(majorC_DGS)
wide_majorC_DGS <- data.table::dcast(setDT(tmpDT), 
                                     name + ct + dx ~ region,
                                     value.var=c('Log2FC','FDR'))

write.csv(wide_majorC_DGS, file.path(projdir, "dxDGS_majorCs_fdr0.1.csv"))


#subC
subC_differentialGene <- readRDS( file.path(projdir, paste("df_dxDEG_fdr0.1_subCs.rds") ) )
subC_dGS_region <- readRDS( file.path(projdir, paste("df_dxDEG_fdr0.1_subCs_region.rds") ) )
cols <- c("name","Log2FC","FDR","groupName","region","subC","dx")

subC_differentialGene$region <- rep("none",dim(subC_differentialGene)[1])

subC_dGS_region$region <- subC_dGS_region$subC
subC_dGS_region$subC <-  gsub(".*\\|(.*?)\\|.*","\\1", subC_dGS_region$groupName)


head(subC_differentialGene)
head(subC_dGS_region)

cols <- colnames(subC_differentialGene)
subC_DGS <- rbind(subC_differentialGene[,cols], subC_dGS_region[,cols])


tmpDT <- as.data.frame(subC_DGS)
wide_subC_DGS <- data.table::dcast(setDT(tmpDT), 
                                   name + subC + dx ~ region,
                                   value.var=c('Log2FC','FDR'))


write.csv(wide_subC_DGS, file.path(projdir, "dxDGS_subCs_fdr0.1.csv"))




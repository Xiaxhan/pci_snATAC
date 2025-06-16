library(ArchR)
projdir <- "projrmSubset_subPeak"
proj <-loadArchRProject(projdir)
dropSubCs <- c("undefined","ast.C6","mg.C8","mg.C10","mg.C15",
               "mg.C1", "mg.C5", "neu.C3", "neu.C4", "odc.C2", "neu.C2", "mg.C2")
proj$subClusters[proj$subClusters %like% "neuron"] <- gsub("neuron","neu",
                                                           proj$subClusters[proj$subClusters %like% "neuron"])
projclean <- proj[!proj$subClusters %in% dropSubCs,]

subCs <- sort(unique(projclean$subClusters))
rm(projclean)
rm(proj)

interaction_grouping <- read.delim("CellPhoneDB/interactions_groups.txt")
colnames(interaction_grouping) <- c("interaction","role")

subCs <- gsub("neu","neuron",subCs)

IN <- subCs[subCs %in% c("neuron.C6","neuron.C7","neuron.C8","neuron.C9")]
EX <- subCs[subCs %like% "neuron"]
EX <- EX[!EX %in% IN] 

neurons <- c(paste("EX", EX, sep = "."), paste("IN", IN, sep = "."))

celltypes <- list() # glia
celltypes$ast <- subCs[subCs %like% "ast"]
celltypes$mg <- subCs[subCs %like% "mg"]
celltypes$odc <- subCs[subCs %like% "odc"]
celltypes$opc <- subCs[subCs %like% "opc"]
celltypes$IN <- IN
celltypes$EX <- EX

ntype <- "subC"
dropSubCs <- gsub("neu","neuron",dropSubCs)
dropSubCs

for(dx in c("AD","bvFTD","PSP_S","Control")){
  print(dx)
  file_pvalue <- paste("degs_analysis_relevant_interactions_result_DEG", ntype, dx, "txt",sep = ".")
  file_sigmeans <- paste("degs_analysis_significant_means_DEG", ntype, dx, "txt",sep = ".")
  file_means <- paste("degs_analysis_means_result_DEG", ntype, dx, "txt",sep = ".")
  flle_decon <- paste("degs_analysis_deconvoluted_result_DEG", ntype, dx,"txt",sep = ".")
  pval.t <- read.delim(file.path(outp, file_pvalue), check.names = FALSE)
  mean.t <- read.delim(file.path(outp, file_means), check.names = FALSE)
  sigmean.t <- read.delim(file.path(outp, file_sigmeans), check.names = FALSE)
  decon.t <- read.delim(file.path(outp, flle_decon), check.names = FALSE)
  
  
  #rm CCI with dropSubCs
  colnames(pval.t)[colnames(pval.t) %like% "neuron.C2"]
  for(rmSubC in dropSubCs){
    dropPairs <- c(colnames(pval.t)[colnames(pval.t) %like% paste(rmSubC,"$",sep="")],
                   colnames(pval.t)[colnames(pval.t) %like% paste(rmSubC,"\\|",sep="")])
    pval.t <- pval.t[,!colnames(pval.t) %in% dropPairs]
    mean.t <- mean.t[,!colnames(mean.t) %in% dropPairs]
    sigmean.t <- sigmean.t[,!colnames(sigmean.t) %in% dropPairs ]
    decon.t <- decon.t[,!colnames(decon.t) %in% rmSubC]
  }
  
  for(iin in IN){
    newNeu <- paste("IN", iin, sep = ".")
    colnames(pval.t)[colnames(pval.t) %like% iin] <- gsub(iin, newNeu, colnames(pval.t)[colnames(pval.t) %like% iin]  )
    colnames(mean.t)[colnames(mean.t) %like% iin] <- gsub(iin,newNeu, colnames(mean.t)[colnames(mean.t) %like% iin]  )
    colnames(sigmean.t)[colnames(sigmean.t) %like% iin] <-  gsub(iin,newNeu, colnames(sigmean.t)[colnames(sigmean.t) %like% iin]  )
    colnames(decon.t)[colnames(decon.t) %like% iin] <- gsub(iin,newNeu, colnames(decon.t)[colnames(decon.t) %like% iin]  )
  }
  colnames(pval.t)[colnames(pval.t) %like% "IN"]
  
  
  colnames(pval.t)[colnames(pval.t) %like% "EX"]
  colnames(pval.t)[colnames(pval.t) %like% "neuron.C2"]
  # rm columns with drop clusters!
  for(ex in c("neuron.C1","neuron.C5")){
    newNeu <- paste("EX",ex, sep = ".")
    colnames(pval.t)[colnames(pval.t) %like% ex] <- gsub(ex,newNeu, colnames(pval.t)[colnames(pval.t) %like% ex]  )
    colnames(mean.t)[colnames(mean.t) %like% ex] <- gsub(ex,newNeu, colnames(mean.t)[colnames(mean.t) %like% ex]  )
    colnames(sigmean.t)[colnames(sigmean.t) %like% ex] <- gsub(ex,newNeu, colnames(sigmean.t)[colnames(sigmean.t) %like% ex]  )
    colnames(decon.t)[colnames(decon.t) %like% ex] <- gsub(ex,newNeu, colnames(decon.t)[colnames(decon.t) %like% ex]  )
  }
  colnames(pval.t)[colnames(pval.t) %like% "EX"]

  if(dx %in% "AD"){ dxName <- "AD vs Control" }
  if(dx %in% "bvFTD"){ dxName <- "bvFTD vs Control" }
  if(dx %in% "PSP_S"){ dxName <- "PSP_S vs Control" }
  if(dx %in% "Control.AD"){ dxName <- "Control vs AD" }
  if(dx %in% "Control.bvFTD"){ dxName <- "Control vs bvFTD" }
  if(dx %in% "Control.PSP_S"){ dxName <- "Control vs PSP_S" }
  if(dx %in% "Control"){ dxName <- "Control vs 3Dx" }
  
  pvalsL[[dxName]] <- pval.t
  meansL[[dxName]] <- mean.t
  sigmeansL[[dxName]] <- sigmean.t
  deconL[[dxName]] <- decon.t
}

newDx <- c("AD vs Control","bvFTD vs Control","PSP_S vs Control", "Control vs 3Dx")

# count nb of cci of each mg subtype
subCsNew <- c(subCs[!subCs %like% "neuron"],neurons)
counts_cci <- NULL
for(dx in c("AD vs Control","bvFTD vs Control","PSP_S vs Control", "Control vs 3Dx")){
  pvals <- pvalsL[[dx]] 
  all_cell_pairs <- colnames(pvals)[12:dim(pvals)[2]]
  for (subC in subCsNew){
    # count CCIs involve in a subC
    CCI_of_subC <- all_cell_pairs[all_cell_pairs %like% subC]
    sigCounts <- rowSums(pvals[, CCI_of_subC])
    length(sigCounts[sigCounts>0])
    dt1 <- data.frame( ct=subC, 
                       sigCCI=length(sigCounts[sigCounts>0]), 
                       group = "cci_involve_subC", 
                       dx =dx)
    counts_cci <- rbind(counts_cci, dt1)
    
    # count CCI of subC with specific cell type 
    for (CT in c("ast","mg","EX","IN","odc","opc")){
      if(!subC %like% CT){
        CCI_pair <- all_cell_pairs[all_cell_pairs %like% subC & all_cell_pairs %like% CT]
      }else{
        otherCTs <- c("ast","mg","EX","IN","odc","opc")
        otherCTs <- otherCTs[!otherCTs %in% CT]
        CCI_pair <- all_cell_pairs[all_cell_pairs %like% subC & 
                                     !all_cell_pairs %like% otherCTs[1] &
                                     !all_cell_pairs %like% otherCTs[2] &
                                     !all_cell_pairs %like% otherCTs[3] &
                                     !all_cell_pairs %like% otherCTs[4] &
                                     !all_cell_pairs %like% otherCTs[5] ]
      }
      sigCounts2 <- rowSums(pvals[, CCI_pair])
      
      dt2 <- data.frame( ct=subC, 
                         sigCCI=length(sigCounts2[sigCounts2>0]), 
                         group = paste("pair_with_", CT, sep = ""), 
                         dx =dx)
      counts_cci <- rbind(counts_cci, dt2)
    }
  }
}
counts_cci[counts_cci$group %like% "EX",]


# barplot
pbar <- list()
counts_cci$dx <- factor(counts_cci$dx, 
                        levels=c("Control vs 3Dx","AD vs Control","bvFTD vs Control", "PSP_S vs Control"))
colors.dx <- ArchR::ArchRPalettes$circus[1:4]
names(colors.dx) <- c("Control vs 3Dx","AD vs Control","bvFTD vs Control", "PSP_S vs Control")

dt1 <- counts_cci[counts_cci$group %in% "cci_involve_subC",]
pbar[["all"]] <- ggplot(dt1, aes(x=ct, y=sigCCI, fill=dx)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.8)+
  #geom_text(aes(label=nb), hjust = 1, size = 2) +
  scale_fill_manual(values = colors.dx)  +
  labs(x="",y="",fill="",title="CCIs involved in ATAC subClusters")+ 
  #coord_flip() +   
  theme(plot.title = element_text(size=10),
        axis.title.x=element_text(size=10,face = "plain"),
        axis.title.y=element_text(size=10,vjust = 2,hjust = 0.5,face = "plain"),
        axis.text.x=element_text(size=10, face="plain", angle = 45, vjust = 1, hjust=1),
        axis.text.y=element_text(size=10,face="plain"),
        legend.position = "right",
        #legend.position = c(1,1),
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        panel.grid =element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=1, colour = "black"))
for (ct in c("ast","mg","EX","IN","odc","opc")){
  dtss <- counts_cci[counts_cci$group %like% ct,]
  pbar[[ct]] <- ggplot(dtss, aes(x=ct, y=sigCCI, fill=dx)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.8)+
    #geom_text(aes(label=nb), hjust = 1, size = 2) +
    scale_fill_manual(values = colors.dx)  +
    labs(x="",y="",fill="",
         title=paste("CCIs with ", ct, " in ATAC subClusters", sep = ""))+ 
    #coord_flip() +   
    theme(plot.title = element_text(size=10),
          axis.title.x=element_text(size=10,face = "plain"),
          axis.title.y=element_text(size=10,vjust = 2,hjust = 0.5,face = "plain"),
          axis.text.x=element_text(size=10, face="plain", angle = 45, vjust = 1, hjust=1),
          axis.text.y=element_text(size=10,face="plain"),
          legend.position = "right",
          #legend.position = c(1,1),
          legend.title=element_blank(),
          legend.text=element_text(size=8),
          panel.grid =element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size=1, colour = "black"))
}
pdf(file.path(outp,"pBar.CCInb_involve_SubC.pdf"), width = 20, height = 15)
ggpubr::ggarrange(plotlist = pbar, nrow=7, ncol=1)
dev.off()


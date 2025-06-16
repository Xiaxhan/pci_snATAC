#Get barcodes of subclusters

library(ArchR)
outp <- "bam"

# load proj
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
unique(projclean$Sample)
projclean 

#---
astC1_cells <- projclean@cellColData[projclean@cellColData$subClusters %in% "ast.C1",] %>% rownames(.)
length(astC1_cells)
write.table(data.frame(astC1_cells),file.path(outp, "barcodes", "astC1_barcodes.txt"), row.names = F, col.names = F, quote = F)

c4_cells <- projclean@cellColData[projclean@cellColData$subClusters %in% "mg.C4",] %>% rownames(.)
length(c4_cells)
write.table(data.frame(c4_cells),file.path(outp, "barcodes", "mgC4_barcodes.txt"), row.names = F, col.names = F, quote = F)

#extract subC barcodes per sample
samples <- projclean@cellColData$Sample %>% unique()
for(subC in c("ast.C1","mg.C4")){
  for (sample in samples){
    subC_per_sample <- projclean@cellColData[projclean@cellColData$subClusters %in% subC &
                                               projclean@cellColData$Sample %in% sample,] %>% rownames(.)
    barcodes <- gsub(sample,"",subC_per_sample) %>% gsub("#","",.)
    barcodes <- paste("CB:Z:",barcodes,sep="")
    
    subCname <- gsub("\\.","",subC)
    fout <- paste(subCname,sample,"barcodes.txt",sep = ".")
    write.table(data.frame(barcodes),file.path(outp, "barcodes", fout), row.names = F, col.names = F, quote = F)
 
    print(paste(subC, sample, "#barcodes=",length(barcodes), sep = " "))
  }
}

fragments_gr <- getFragmentsFromProject(projclean, cellNames = subcluster_cells)

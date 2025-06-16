library(Signac)
library(Seurat)
library(SeuratDisk)
library(dplyr)
outp <- "CellPhoneDB/output"

# sobj of dowmsampled cellset 
sobj.ct <- LoadH5Seurat(file.path(outp, "tmp.Control.scATAC_clean.h5seurat"))
sobj.ad <- LoadH5Seurat(file.path(outp, "tmp.AD.scATAC_clean.h5seurat"))
sobj.ftd <- LoadH5Seurat(file.path(outp, "tmp.bvFTD.scATAC_clean.h5seurat"))
sobj.psp <- LoadH5Seurat(file.path(outp, "tmp.PSP_S.scATAC_clean.h5seurat"))


sobj.ct_ad <- merge(x=sobj.ct, y=sobj.ad,  add.cell.ids = NULL)
sobj.ct_ftd <- merge(x=sobj.ct, y=sobj.ftd,  add.cell.ids = NULL)
sobj.ct_psp <- merge(x=sobj.ct, y=sobj.psp,  add.cell.ids = NULL)
rm(sobj.ct)
rm(sobj.ad)
rm(sobj.ftd)
rm(sobj.psp)

sobjlist<-list("AD"=sobj.ct_ad, "bvFTD"=sobj.ct_ftd, "PSP_S"=sobj.ct_psp)

#----Dx DEG in cell types
dtDEG.dx2ct <- NULL
for (dx in c("AD","bvFTD","PSP_S")){
  print(dx)
  sobj <- sobjlist[[dx]]
  Idents(object = sobj) <- sobj$majorClusters
  table(Idents(sobj))
  for (celltype in unique(sobj$majorClusters)){
    print(celltype)
    dxDEG <- FindMarkers(sobj, ident.1 = dx, group.by = 'Clinical.Dx', subset.ident = celltype)
    n_rows <- dim(dxDEG)[1]
    dxDEG <- cbind(dxDEG, dx=rep(dx, n_rows), celltype = rep(celltype, n_rows))
    dtDEG.dx2ct <- rbind(dtDEG.dx2ct, dxDEG)
  }
}

for (dx in c("AD","bvFTD","PSP_S")){
  dx.DEG <- dtDEG.dx2ct[dtDEG.dx2ct$dx %in% dx,]
  colnames(dx.DEG) <- c("P.value","logFC","pct.1","pct.2","adj.P.Val","dx","clusters")
  dx.DEG <- cbind(dx.DEG, gene=rownames(dx.DEG))
  up.DEG <- dx.DEG[dx.DEG$logFC>0,c("clusters","gene")]
  fname <- paste("DEG",dx,"downsample.txt",sep = ".")
  write.table(up.DEG, file.path(outp, fname), quote = F, sep = "\t", row.names = F)
}



#----------Dx DEG in subtypes
subClusters <- sort(unique(sobjlist$AD$subClusters))
dtDEG.dx2ct <- NULL
for (dx in c("AD","bvFTD","PSP_S")){
  #print(dx)
  sobj <- sobjlist[[dx]]
  Idents(object = sobj) <- sobj$subClusters
  table(Idents(sobj))
  for (celltype in subClusters){
    print(paste(dx,celltype))
    dxDEG <- FindMarkers(sobj, ident.1 = dx, group.by = 'Clinical.Dx', subset.ident = celltype)
    n_rows <- dim(dxDEG)[1]
    dxDEG <- cbind(dx=rep(dx, n_rows), celltype = rep(celltype, n_rows), gene = rownames(dxDEG), dxDEG)
    rownames(dxDEG) <- NULL
    dtDEG.dx2ct <- rbind(dtDEG.dx2ct, dxDEG)
  }
}
write.csv(dtDEG.dx2ct, file.path(outp, "alldeg.subC.3dx.csv")) # raw WO filtering

dtDEG <- dtDEG.dx2ct

dtDEG.dx2ct <- read.csv(file.path(outp, "alldeg.subC.3dx.csv"), header = T)
#dx up
log2fc_cutoff <- 0.25  
min(abs(dtDEG.dx2ct$avg_log2FC))
#[1] 0.2500004

for (dx in c("AD","bvFTD","PSP_S")){
  dx.DEG <- dtDEG.dx2ct[dtDEG.dx2ct$dx %in% dx,]
  up.DEG <- dx.DEG[dx.DEG$avg_log2FC>log2fc_cutoff & dx.DEG$p_val < 0.05 & dx.DEG$p_val_adj <0.1, 
                   c("celltype","gene")]
  fname <- paste("DEG.subC",dx,"downsample.txt",sep = ".")
  write.table(up.DEG, file.path(outp, fname), quote = F, sep = "\t", row.names = F)
  #save for reading
  tmpDEG <- dx.DEG[dx.DEG$avg_log2FC>log2fc_cutoff & dx.DEG$p_val < 0.05 & dx.DEG$p_val_adj <0.1,]
  fname <- paste("checkDEG.subC",dx,"downsample.txt",sep = ".")
  write.table(tmpDEG, file.path(outp, fname), quote = F, sep = "\t", row.names = F)
}

#ct up
for (dx in c("AD","bvFTD","PSP_S")){
  dx.DEG <- dtDEG.dx2ct[dtDEG.dx2ct$dx %in% dx,]
  up.DEG <- dx.DEG[dx.DEG$avg_log2FC<log2fc_cutoff & dx.DEG$p_val < 0.05 & dx.DEG$p_val_adj <0.1, 
                   c("celltype","gene")]
  fname <- paste("DEG.subC.Control",dx,"downsample.txt",sep = ".")
  write.table(up.DEG, file.path(outp, fname), quote = F, sep = "\t", row.names = F)
  #save for reading
  tmpDEG <- dx.DEG[dx.DEG$avg_log2FC<log2fc_cutoff & dx.DEG$p_val < 0.05 & dx.DEG$p_val_adj <0.1,]
  fname <- paste("checkDEG.subC.Control",dx,"downsample.txt",sep = ".")
  write.table(tmpDEG, file.path(outp, fname), quote = F, sep = "\t", row.names = F)
}

#----------Control DEG in subcs!!!! compare with all 3dx
library(Signac)
library(Seurat)
library(SeuratDisk)
library(dplyr)
outp <- "CellPhoneDB/output"

# sobj of dowmsampled cellset 
sobj.ct <- LoadH5Seurat(file.path(outp, "tmp.Control.scATAC_clean.h5seurat"))
sobj.ad <- LoadH5Seurat(file.path(outp, "tmp.AD.scATAC_clean.h5seurat"))
sobj.ftd <- LoadH5Seurat(file.path(outp, "tmp.bvFTD.scATAC_clean.h5seurat"))
sobj.psp <- LoadH5Seurat(file.path(outp, "tmp.PSP_S.scATAC_clean.h5seurat"))

sobjlist<-list("Control"=sobj.ct,"AD"=sobj.ad, "bvFTD"=sobj.ftd, "PSP_S"=sobj.psp)

subClusters <- sort(unique(sobjlist$AD$subClusters))


dtDEG.ct <- NULL
for (celltype in subClusters){
  print(celltype)
  ssct <- subset(x = sobj.ct, subset = subClusters == celltype)
  merge_sobj <- ssct
  print("merging sobj subsets...")
  for (dx in c("AD","bvFTD","PSP_S")){
    print(dx)
    ssdx <- subset(x = sobjlist[[dx]], subset = subClusters == celltype)
    merge_sobj <- merge(x=merge_sobj, y=ssdx,  add.cell.ids = NULL)
  }
 
  dxDEG <- FindMarkers(merge_sobj, ident.1 = "Control", group.by = 'Clinical.Dx')
  n_rows <- dim(dxDEG)[1]
  dxDEG <- cbind(dx=rep("Control", n_rows), celltype = rep(celltype, n_rows), gene = rownames(dxDEG), dxDEG)
  rownames(dxDEG) <- NULL
  dtDEG.ct <- rbind(dtDEG.ct, dxDEG)
}
write.csv(dtDEG.ct, file.path(outp, "alldeg.subC.control_vs_3dx.csv")) # raw WO filtering


log2fc_cutoff <- 0.25 # 0 = 0.25
up.DEG <- dtDEG.ct[dtDEG.ct$avg_log2FC>log2fc_cutoff & dtDEG.ct$p_val < 0.05 & dtDEG.ct$p_val_adj <0.1, 
                 c("celltype","gene")]
fname <-"DEG.subC.Control.downsample.txt"
write.table(up.DEG, file.path(outp, fname), quote = F, sep = "\t", row.names = F)
#save for reading
tmpDEG <- dtDEG.ct[dtDEG.ct$avg_log2FC>log2fc_cutoff & dtDEG.ct$p_val < 0.05 & dtDEG.ct$p_val_adj <0.1,]
fname <- "checkDEG.subC.Control.downsample.txt"
write.table(tmpDEG, file.path(outp, fname), quote = F, sep = "\t", row.names = F)


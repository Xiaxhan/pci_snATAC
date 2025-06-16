## Validate the identified Enhancers

library(ggplot2)
library(ArchR)
addArchRThreads(threads = 10)

projdir <- "projrmSubset_subPeak" 
dt.Enhancer.majorGroupExist <- readRDS(file.path(getwd(), "outdir_projrmSubset_subPeak/CRE_projrmSubset_subPeak","subC.enhancer_majorGroupExist.rds"))
#!!!!!!!!!!!

dt.exist <- dt.Enhancer.majorGroupExist
dt.exist[dt.exist > 0] <- 1
dt.exist <- as.data.frame(dt.exist)
dt.exist <- cbind(dt.exist, rowsum = rowSums(dt.exist))
dt1 <- as.data.frame(table(dt.exist$rowsum))

#save enh by cell types - all enh - unique enh
dt.exist <- dt.Enhancer.majorGroupExist
dt.exist[dt.exist > 0] <- 1
dt.exist <- cbind(dt.exist, rowsum = rowSums(dt.exist))

as.data.frame(table(dt.exist[,"rowsum"]))

dt.uniqueEnh <- dt.exist[dt.exist[,"rowsum"]==1,]
dt.commonAllEnh <- dt.exist[dt.exist[,"rowsum"]==5,]

sum.Enhancer.majorGroupExist <- cbind(dt.Enhancer.majorGroupExist, 
                                     rowsum = rowSums(dt.Enhancer.majorGroupExist))
sum.Enhancer.majorGroupExist[rownames(dt.uniqueEnh),]

# control
dt.labels_dx <- dt.Enhancer.majorGroupExist
for (j in 1:5){
  cellcol <- dt.Enhancer.majorGroupExist[,j]
  flag_control<-as.numeric(grepl("1...$",cellcol)) # in control 1xxx
  flag_AD <- as.numeric(grepl("1..$",cellcol)) # in AD x1xx
  flag_FTD <- as.numeric(grepl("1.$",cellcol))  #xx1x
  flag_PSP <- as.numeric(grepl("1$",cellcol))  #xx1x
  
  dt.labels_dx <- cbind(dt.labels_dx, flag_control, flag_AD, flag_FTD, flag_PSP)
}
colnames(dt.labels_dx) <- c(colnames(dt.labels_dx)[1:5], paste("a", c("C","A","F","P"),sep=""),
                            paste("m", c("C","A","F","P"),sep=""), paste("n", c("C","A","F","P"),sep=""),
                            paste("d", c("C","A","F","P"),sep=""), paste("p", c("C","A","F","P"),sep=""))
head(dt.labels_dx)

controlCols <- c("aC",'mC',"nC","dC","pC")
dt.labels_dx <- cbind(dt.labels_dx, sum_control = rowSums(dt.labels_dx[,controlCols]), row_sum = rowSums(dt.labels_dx[,6:25]) )
head(dt.labels_dx)

# if focus on Enhancers which detect in control (no matter whether it exist in Dx)
dt.control <- dt.labels_dx[,controlCols]
dt.table <- as.data.frame(table(dt.labels_dx[,"sum_control"]))[-1,]
sum(dt.table$Freq)
dt.table <- cbind(dt.table, freq=dt.table$Freq/sum(dt.table$Freq))
dt.table


# drop E(P-E) not in subtype-diagnose specific group
dim(dt.Enhancer.majorGroupExist)  #223710
filterdE <- rownames(dt.exist[dt.exist[,"rowsum"]==0,])
dt.newEnhancer <- dt.Enhancer.majorGroupExist[rownames(dt.Enhancer.majorGroupExist) %ni% filterdE,]

gr.enhancer.celltype <-list()
gr.enhancer.celltypeUnique <-list()
for(cell in colnames(dt.newEnhancer)){
  #all enh
  dt <- dt.newEnhancer[,cell]
  enhancers <- names(dt[dt!=0])
  length(enhancers)
  cEnhan <- unlist(strsplit(enhancers,"-"))
  idx_chrR <- seq(1, length(cEnhan)-2, by = 3)
  idx_startR <- seq(2, length(cEnhan)-1, by = 3)
  idx_endR <- seq(3, length(cEnhan), by = 3)
  chrR <- cEnhan[idx_chrR]
  start.R <- as.numeric(cEnhan[idx_startR])
  end.R <- as.numeric(cEnhan[idx_endR])
  
  gr.enhancer.celltype[[cell]] <- GRanges(
    seqnames = Rle(chrR),
    ranges = IRanges(start.R, end = end.R))
  
  #unique enh
  dt <- dt.uniqueEnh[,cell]
  enhancers <- names(dt[dt!=0])
  length(enhancers)
  cEnhan <- unlist(strsplit(enhancers,"-"))
  idx_chrR <- seq(1, length(cEnhan)-2, by = 3)
  idx_startR <- seq(2, length(cEnhan)-1, by = 3)
  idx_endR <- seq(3, length(cEnhan), by = 3)
  chrR <- cEnhan[idx_chrR]
  start.R <- as.numeric(cEnhan[idx_startR])
  end.R <- as.numeric(cEnhan[idx_endR])
  
  gr.enhancer.celltypeUnique[[cell]] <- GRanges(
    seqnames = Rle(chrR),
    ranges = IRanges(start.R, end = end.R))
}

cEnhan <- rownames(dt.newEnhancer)
cEnhan <- unlist(strsplit(cEnhan,"-"))
idx_chrR <- seq(1, length(cEnhan)-2, by = 3)
idx_startR <- seq(2, length(cEnhan)-1, by = 3)
idx_endR <- seq(3, length(cEnhan), by = 3)
chrR <- cEnhan[idx_chrR]
start.R <- as.numeric(cEnhan[idx_startR])
end.R <- as.numeric(cEnhan[idx_endR])

gr.enhancer <- GRanges(
  seqnames = Rle(chrR),
  ranges = IRanges(start.R, end = end.R))

gr.enhancer

#--------------------------overlap results
refEnh <- list()
overlapList<-list()
overlap.celltype<-list()

##---------------------------------------------annotate 1) NottScience
path_epiAnno <- "/home/xhan/data_cellrangerOut/epiResource/NottScience"
cell_epiAnno <- c("LHX2", "NeuN", "Olig2", "PU1") #ast neu oligo mg
file_epiAnno <- paste(cell_epiAnno, "_enhancers.sorted.bed", sep = "")
dt.epiAnno <- NULL
for (file in file_epiAnno){
  print(file)
  dt <- read.csv(file.path(path_epiAnno,file), sep = "\t", header = F)
  dt <- cbind(dt, celltype = rep(strsplit(file, "_")[[1]][1], dim(dt)[1]))
  head(dt)
  dt.epiAnno <- rbind(dt.epiAnno, dt)
}
colnames(dt.epiAnno) <- c("chr","start","end","celltype")
dim(dt.epiAnno)
head(dt.epiAnno)

gr.epiAnno <- GRanges(
  seqnames = Rle(dt.epiAnno$chr),
  ranges = IRanges(dt.epiAnno$start, end = dt.epiAnno$end),
  celltype = dt.epiAnno$celltype)
gr.epiAnno

gr.epiAnno.unique <- unique(gr.epiAnno)
gr.epiAnno.unique$hits <- countOverlaps(gr.epiAnno.unique, gr.epiAnno, type="equal")

gr.epiAnno.common <- gr.epiAnno.unique[gr.epiAnno.unique$hits!=1]


#check overlapping
o <- findOverlaps(unique(gr.epiAnno), gr.enhancer)
nb_annoEnhan <- length(unique(gr.epiAnno))
nb_annoEnhan #316075
nb_myEnhan <- length(gr.enhancer)
nb_myEnhan

gr.overlap_in_ourEn <- unique(gr.enhancer[subjectHits(o)])

overlap_in_ourEn <- unique(gr.enhancer[subjectHits(o)])
nb_overlap <- length(overlap_in_ourEn)
nb_overlap

nb_overlap / nb_myEnhan #0.2071632

overlap_in_ourEn <- paste(seqnames(overlap_in_ourEn), overlap_in_ourEn@ranges, sep = "-")

overlapList[["NottScience"]] <- overlap_in_ourEn 
refEnh[["NottScience"]] <- nb_annoEnhan

## overlap analysis - split by cell types - use enhancers exist in cell type (no matter whether it exist in other cell types)
dt <- c("celltype","nb_ours","nb_anno","overlap_ours","overlap_anno","freq.overlap_ours","freq.overlap_anno")
names(cell_epiAnno) <- c("ast", "neuron", "odc", "mg")
for (cell in names(gr.enhancer.celltype)){
  #cell  <- "ast"
  if(cell %in% names(cell_epiAnno)){
    annoteType <- cell_epiAnno[cell]
    
    gr.ours <- gr.enhancer.celltype[[cell]]
    gr.anno <- gr.epiAnno[gr.epiAnno$celltype %in% annoteType]
    
    o <- findOverlaps(gr.ours, gr.anno)
    o
    
    overlap.ours <- unique(gr.ours[queryHits(o)])
    overlap.anno <- unique(gr.anno[subjectHits(o)])
    
    freq.overlap.ours <- length(overlap.ours)/ length(gr.ours)
    freq.overlap.anno <- length(overlap.anno)/ length(gr.anno)
      
      
    saveLine <- c(cell, length(gr.ours), length(gr.anno), length(overlap.ours), length(overlap.anno), freq.overlap.ours, freq.overlap.anno)
    
    dt<-rbind(dt, saveLine)
  }
}
colnames(dt) <- dt[1,]
rownames(dt) <- dt[,1]
dt<-dt[-1,-1]
dt 
overlap.celltype[["withNottScience_exactOverlap"]] <- dt

#among these overlapped Enh(overlap_in_ourEn, 45627), how they distri among our cell types?
#aaaa<-gr.overlap_in_ourEn
#gr.overlap_in_ourEn <- aaaa
for (cell in names(gr.enhancer.celltype)){
    gr.cell <- gr.enhancer.celltype[[cell]]
    m <- findMatches(gr.overlap_in_ourEn, gr.cell)
    cellExist <- rep(0, length(gr.overlap_in_ourEn))
    cellExist[queryHits(m)] <- 1
    
    values(gr.overlap_in_ourEn) <- cbind(values(gr.overlap_in_ourEn), DataFrame(cellExist))  
}
colnames(values(gr.overlap_in_ourEn)) <- names(gr.enhancer.celltype)
gr.overlap_in_ourEn

out <- NULL
toCount <- data.frame(values(gr.overlap_in_ourEn)@listData)
toCount <- cbind(toCount, sum=rowSums(toCount))
unique <- toCount[toCount$sum==1,]
for (cell in names(gr.enhancer.celltype)){
  dt <- values(gr.overlap_in_ourEn)[,cell]
  freq <- length(dt[dt!=0])/ length(gr.overlap_in_ourEn)
  out <- cbind(out, c(cell, length(dt[dt!=0]), freq))
}
out2 <- NULL
for (cell in names(gr.enhancer.celltype)){
  ## count those celltype unique
  dt2 <- unique[,cell]
  freq2 <- length(dt2[dt2!=0])/ length(gr.overlap_in_ourEn)
  out2 <- cbind(out2, c(paste(cell,"Unique", sep = ""), 
                      length(dt2[dt2!=0]), freq2))
}
out <- rbind(out, out2)
#colnames(out) <- out[1,]
rownames(out) <- rep(c("celltype", "nb_enh", "freq"),2)
out

overlap.celltype[["withNottScience_DistriAmongOurCell"]] <- out

## overlap analysis - split by cell types - only cell type unique enhancers (also for annotate)
dt <- c("celltype","nb_ours","nb_anno","overlap_ours","overlap_anno","freq.overlap_ours","freq.overlap_anno")
names(cell_epiAnno) <- c("ast", "neuron", "odc", "mg")
for (cell in names(gr.enhancer.celltype)){
  #cell  <- "ast"
  if(cell %in% names(cell_epiAnno)){
    annoteType <- cell_epiAnno[cell]
    
    gr.ours <- gr.enhancer.celltypeUnique[[cell]]
    gr.anno <- gr.epiAnno[gr.epiAnno$celltype %in% annoteType]
    gr.anno <- setdiff(gr.anno, gr.epiAnno.common) # drop the celltype common annotatedEnhs
    
    o <- findOverlaps(gr.ours, gr.anno)
    o
    
    overlap.ours <- unique(gr.ours[queryHits(o)])
    overlap.anno <- unique(gr.anno[subjectHits(o)])
    
    freq.overlap.ours <- length(overlap.ours)/ length(gr.ours)
    freq.overlap.anno <- length(overlap.anno)/ length(gr.anno)
    
    
    saveLine <- c(cell, length(gr.ours), length(gr.anno), length(overlap.ours), length(overlap.anno), freq.overlap.ours, freq.overlap.anno)
    
    dt<-rbind(dt, saveLine)
  }
}
colnames(dt) <- dt[1,]
rownames(dt) <- dt[,1]
dt<-dt[-1,-1]
dt
overlap.celltype[["withNottScience_exactOverlap_cellUniqueEnh"]] <- dt

##---------------------------------------------annotate 2) ENCODE
## cCREs
path_encode <- "/home/xhan/data_cellrangerOut/epiResource/ENCODE"
file_encode <- c("ENCFF198KYT_bipolarNeuron_cCREs.bed", "ENCFF815WRM_ast_cCREs.bed") # X EXNCFF465ARD_frontalCortex_cCREs.bed - only "Low-DNase"  "DNase-only"????
for (file in file_encode){
  table.encode <- read.csv(file.path(path_encode, file), sep = "\t", header = F)
  head(table.encode)
  unique(table.encode$V10)
  unique(table.encode$V11)
  unique(table.encode$V9)
  
  enhancer_label.encode <- c("pELS", "dELS", "dELS,CTCF-bound", "pELS,CTCF-bound")
  enhancer.encode <- table.encode[table.encode$V10 %in% enhancer_label.encode,]
  dim(enhancer.encode)[1] #74853
  tail(enhancer.encode)
  gr.enhancer.encode <- GRanges(
    seqnames = Rle(enhancer.encode$V1),
    ranges = IRanges(enhancer.encode$V2, end = enhancer.encode$V3),
    type = enhancer.encode$V10)
  
  gr.enhancer.encode
  unique(gr.enhancer.encode)
  #check overlapping
  o2 <- findOverlaps(gr.enhancer.encode, gr.enhancer)
  nb_annoEnhan <- length(unique(gr.enhancer.encode))
  nb_annoEnhan #74853
  overlap_in_ourEn <- unique(gr.enhancer[subjectHits(o2)])
  nb_overlap <- length(overlap_in_ourEn)
  nb_overlap
  
  nb_overlap / 218414 #0.06767881
  
  overlap_in_ourEn <- paste(seqnames(overlap_in_ourEn), overlap_in_ourEn@ranges, sep = "-")
  length(overlap_in_ourEn)

  name <- strsplit(file, "_")[[1]][2]
  name <- paste("ENCODE",name, "cCREs", sep = "_")
  overlapList[[name]] <- overlap_in_ourEn 
  refEnh[[name]] <- nb_annoEnhan
}

## chromState
file_encode <- c("ENCFF791URB_neuron_chromState.bed", "ENCFF395QLP_ast_chromState.bed")
for (file in file_encode){
  table.encode <- read.csv(file.path(path_encode, file), sep = "\t", header = F)
  dim(table.encode)
  head(table.encode)
  unique(table.encode$V4)
  
  enhancer_label.encode <- c("EnhG1", "EnhG2", "EnhA1", "EnhA2", "EnhWk", "EnhBiv")
  enhancer.encode <- table.encode[table.encode$V4 %in% enhancer_label.encode,]
  dim(enhancer.encode)[1] #145291, 163851
  tail(enhancer.encode)
  gr.enhancer.encode <- GRanges(
    seqnames = Rle(enhancer.encode$V1),
    ranges = IRanges(enhancer.encode$V2, end = enhancer.encode$V3),
    type = enhancer.encode$V4)
  
  gr.enhancer.encode
  unique(gr.enhancer.encode)
  #check overlapping
  o2 <- findOverlaps(gr.enhancer.encode, gr.enhancer)
  nb_annoEnhan <- length(unique(gr.enhancer.encode))
  nb_annoEnhan #145291, 163851
  overlap_in_ourEn <- unique(gr.enhancer[subjectHits(o2)])
  nb_overlap <- length(overlap_in_ourEn)
  nb_overlap # 35574 # 36361
  
  nb_overlap / 218414 #0.1628742, 0.1664774
  
  overlap_in_ourEn <- paste(seqnames(overlap_in_ourEn), overlap_in_ourEn@ranges, sep = "-")
  length(overlap_in_ourEn)
  
  name <- strsplit(file, "_")[[1]][2]
  name <- paste("ENCODE",name, "chromState", sep = "_")
  overlapList[[name]] <- overlap_in_ourEn 
  refEnh[[name]] <- nb_annoEnhan
}

##---------------------------------------------annotate 3) FANTOM5
path_fantom5<-"/home/xhan/data_cellrangerOut/epiResource/FANTOM5"
file_fantom5<-c("enhancer_ast_atleast0.05percent_9666.bed", "enhancer_brain_atleast0.05percent_29106.bed",
               "enhancer_neuron_atleast0.05percent_16806.bed")
for (file in file_fantom5){
  table.encode <- read.csv(file.path(path_fantom5, file), sep = "\t", header = F)
  dim(table.encode)
  head(table.encode)
  
  gr.enhancer.encode <- GRanges(
    seqnames = Rle(table.encode$V1),
    ranges = IRanges(table.encode$V2, end = table.encode$V3))
  
  gr.enhancer.encode
  unique(gr.enhancer.encode)
  #check overlapping
  o2 <- findOverlaps(gr.enhancer.encode, gr.enhancer)
  nb_annoEnhan <- length(unique(gr.enhancer.encode))
  overlap_in_ourEn <- unique(gr.enhancer[subjectHits(o2)])
  nb_overlap <- length(overlap_in_ourEn)
  
  overlap_in_ourEn <- paste(seqnames(overlap_in_ourEn), overlap_in_ourEn@ranges, sep = "-")
  length(overlap_in_ourEn)
  
  name <- strsplit(file, "_")[[1]][2]
  name <- paste("FANTOM",name,sep = "_")
  overlapList[[name]] <- overlap_in_ourEn 
  refEnh[[name]] <- nb_annoEnhan
}

##---------------------------------------------4) ABC
path_ABC <- "/home/xhan/data_cellrangerOut/epiResource/ABC"
file_ABC <- c("EnhancerPredictionsFull_ast.txt", "EnhancerPredictionsFull_biploarNeuron.txt")
for (file in file_ABC){
  table.encode <- read.csv(file.path(path_ABC, file), sep = "\t", header = T)
  dim(table.encode)
  head(table.encode)
  
  gr.enhancer.encode <- GRanges(
    seqnames = Rle(table.encode$chr),
    ranges = IRanges(table.encode$start, end = table.encode$end))
  
  gr.enhancer.encode
  unique(gr.enhancer.encode)
  #check overlapping
  o2 <- findOverlaps(gr.enhancer.encode, gr.enhancer)
  nb_annoEnhan <- length(unique(gr.enhancer.encode))
  overlap_in_ourEn <- unique(gr.enhancer[subjectHits(o2)])
  nb_overlap <- length(overlap_in_ourEn)
  
  overlap_in_ourEn <- paste(seqnames(overlap_in_ourEn), overlap_in_ourEn@ranges, sep = "-")
  length(overlap_in_ourEn)
  
  name <- strsplit(file, "_")[[1]][2]
  name <- paste("ABC",name,sep = "_")
  name
  overlapList[[name]] <- overlap_in_ourEn 
  refEnh[[name]] <- nb_annoEnhan
}

##---------------------------------------------annotate 5) SwarupNC
file<-"/home/xhan/data_cellrangerOut/epiResource/SwarupNC_glCREs.csv" # paper suppleTable4
table.encode <- read.csv(file, sep = "\t", header = T)
dim(table.encode)
head(table.encode)
length(table.encode$cCRE)
length(unique(table.encode$cCRE)) #48880
unique(table.encode$group)


tmp <- unlist(strsplit(unique(table.encode$cCRE),split = "-"))
idx_chr <- seq(1, length(tmp)-2, by = 3)
idx_start <- seq(2, length(tmp)-1, by = 3)
idx_end <- seq(3, length(tmp), by = 3)
chr <- tmp[idx_chr]
start<- as.numeric(tmp[idx_start])
end<- as.numeric(tmp[idx_end])
gr.enhancer.encode <- GRanges(seqnames = Rle(chr), ranges = IRanges(start, end = end))
gr.enhancer.encode

unique(table.encode$group)
unique(table.encode$celltype)

##### nnnew
tmp <- unlist(strsplit(table.encode$cCRE,split = "-"))
idx_chr <- seq(1, length(tmp)-2, by = 3)
idx_start <- seq(2, length(tmp)-1, by = 3)
idx_end <- seq(3, length(tmp), by = 3)
chr <- tmp[idx_chr]
start<- as.numeric(tmp[idx_start])
end<- as.numeric(tmp[idx_end])
gr.epiAnno <- GRanges(seqnames = Rle(chr), ranges = IRanges(start, end = end),
                              celltype = table.encode$celltype, group = table.encode$group)
gr.epiAnno

gr.epiAnno.unique <- unique(gr.epiAnno)
gr.epiAnno.unique$hits <- countOverlaps(gr.epiAnno.unique, gr.epiAnno, type="equal")
gr.epiAnno.common <- gr.epiAnno.unique[gr.epiAnno.unique$hits!=1]


unique(table.encode$celltype)
dt <- data.frame(cCRE=table.encode$cCRE, celltype=table.encode$celltype)
dt <- unique(dt)
dt <- as.data.frame(table(dt$cCRE))
glCRE.rep <- dt[dt$Freq!=1,]
unique(glCRE.rep$Freq)
dim(glCRE.rep) #4680 common celltype-common cCRE

glCRE.allcommon <- glCRE.rep[glCRE.rep$Freq==5,] ##not have 6!!!!
dim(glCRE.allcommon) #6

table.encode[table.encode$cCRE %in% glCRE.allcommon$Var1[1],]


count.dt <- NULL
### Enhancer distri in cell types
for (cell in unique(table.encode$celltype)){
  cell.enh <- table.encode[table.encode$celltype %in% cell,]
  nb_cellEnh <- dim(cell.enh)[1]
  cell.enh.unique <- cell.enh[cell.enh$cCRE  %ni% glCRE.rep$Var1,]
  nb_cellUniqueEnh <- dim(cell.enh.unique)[1]
  count.dt <- rbind(count.dt, c(cell, nb_cellEnh, nb_cellUniqueEnh))
}
count.dt

#check overlapping
o2 <- findOverlaps(gr.enhancer.encode, gr.enhancer)
nb_annoEnhan <- length(unique(gr.enhancer.encode))
gr.overlap_in_ourEn <- unique(gr.enhancer[subjectHits(o2)])
nb_overlap <- length(gr.overlap_in_ourEn)

overlap_in_ourEn <- paste(seqnames(gr.overlap_in_ourEn), gr.overlap_in_ourEn@ranges, sep = "-")
length(overlap_in_ourEn)
overlapList[["SwarupNC"]] <- overlap_in_ourEn 
refEnh[["SwarupNC"]] <- nb_annoEnhan
str(overlapList)

#among these overlapped Enh(overlap_in_ourEn, 19,185), how much of them are from our non-opc cell types?

for (cell in names(gr.enhancer.celltype)){
  gr.cell <- gr.enhancer.celltype[[cell]]
  m <- findMatches(gr.overlap_in_ourEn, gr.cell)
  cellExist <- rep(0, length(gr.overlap_in_ourEn))
  cellExist[queryHits(m)] <- 1
  
  values(gr.overlap_in_ourEn) <- cbind(values(gr.overlap_in_ourEn), DataFrame(cellExist))  
}
colnames(values(gr.overlap_in_ourEn)) <- names(gr.enhancer.celltype)
gr.overlap_in_ourEn

toCount <- data.frame(values(gr.overlap_in_ourEn)@listData)
toCount <- cbind(toCount, sum=rowSums(toCount))
unique <- toCount[toCount$sum==1,]

out <- NULL
for (cell in names(gr.enhancer.celltype)){
  dt <- values(gr.overlap_in_ourEn)[,cell]
  freq <- length(dt[dt!=0])/ length(gr.overlap_in_ourEn)
  out <- cbind(out, c(cell, length(dt[dt!=0]), freq))
}
out2<-NULL
for (cell in names(gr.enhancer.celltype)){
  ## count those celltype unique
  dt2 <- unique[,cell]
  freq2 <- length(dt2[dt2!=0])/ length(gr.overlap_in_ourEn)
  out2 <- cbind(out2, c(paste(cell,"Unique", sep = ""), 
                      length(dt2[dt2!=0]), freq2))
  
}
out <- rbind(out, out2)
rownames(out) <- rep(c("celltype", "nb_enh", "freq"), 2)
out
overlap.celltype[["withSwarupNC_DistriAmongOurCell"]] <- out

## overlap analysis - split by cell types - use enhancers exist in cell type (no matter whether it exist in other cell types)
dt <- c("celltype","nb_ours","nb_anno","overlap_ours","overlap_anno","freq.overlap_ours","freq.overlap_anno")
cell_epiAnno <- unique(gr.epiAnno$celltype)
cell_epiAnno
names(cell_epiAnno) <- c("ast", "neuron", "neuron","mg", "odc","opc")

for (cell in names(gr.enhancer.celltype)){
  if(cell %ni% "neuron"){
    annoteType <- cell_epiAnno[cell]
  }else{
    annoteType <- c("EX","INH")
  }
    gr.ours <- gr.enhancer.celltype[[cell]]
    gr.anno <- gr.epiAnno[gr.epiAnno$celltype %in% annoteType]
    
    o <- findOverlaps(gr.ours, gr.anno)
    o
    
    overlap.ours <- unique(gr.ours[queryHits(o)])
    overlap.anno <- unique(gr.anno[subjectHits(o)])
    
    freq.overlap.ours <- length(overlap.ours)/ length(gr.ours)
    freq.overlap.anno <- length(overlap.anno)/ length(gr.anno)
    
    
    saveLine <- c(cell, length(gr.ours), length(gr.anno), length(overlap.ours), length(overlap.anno), freq.overlap.ours, freq.overlap.anno)
    
    dt<-rbind(dt, saveLine)
  
}
colnames(dt) <- dt[1,]
rownames(dt) <- dt[,1]
dt<-dt[-1,-1]
dt 
overlap.celltype[["withSwarupNC_exactOverlap"]] <- out

## overlap analysis - split by cell types - only cell type unique enhancers (also for annotate)
dt <- c("celltype","nb_ours","nb_anno","overlap_ours","overlap_anno","freq.overlap_ours","freq.overlap_anno")
names(cell_epiAnno) <- c("ast", "neuron", "neuron","mg", "odc","opc")
for (cell in names(gr.enhancer.celltype)){
  #cell  <- "ast"
  if(cell %in% names(cell_epiAnno)){
    annoteType <- cell_epiAnno[cell]
    
    gr.ours <- gr.enhancer.celltypeUnique[[cell]]
    gr.anno <- gr.epiAnno[gr.epiAnno$celltype %in% annoteType]
    gr.anno <- setdiff(gr.anno, gr.epiAnno.common) # drop the celltype common annotatedEnhs
    
    o <- findOverlaps(gr.ours, gr.anno)
    o
    
    overlap.ours <- unique(gr.ours[queryHits(o)])
    overlap.anno <- unique(gr.anno[subjectHits(o)])
    
    freq.overlap.ours <- length(overlap.ours)/ length(gr.ours)
    freq.overlap.anno <- length(overlap.anno)/ length(gr.anno)
    
    
    saveLine <- c(cell, length(gr.ours), length(gr.anno), length(overlap.ours), length(overlap.anno), freq.overlap.ours, freq.overlap.anno)
    
    dt<-rbind(dt, saveLine)
  }
}
colnames(dt) <- dt[1,]
rownames(dt) <- dt[,1]
dt<-dt[-1,-1]
dt
overlap.celltype[["withSwarupNC_exactOverlap_cellUniqueEnh"]] <- dt


## if focus on control-only!
forControl <- dt.Enhancer.majorGroupExist
forControl <- cbind(forControl, sum=rowSums(forControl))
forControl <- forControl[forControl[,"sum"] %% 1000 ==0,] # in all cell types are control

gr.enhancer.control <-list()
gr.enhancer.controlUnique <- list()
control.uniqueEnh <- forControl[forControl[,"sum"]==1000,]
for(cell in colnames(dt.newEnhancer)){
  #all enh
  dt <- forControl[,cell]
  enhancers <- names(dt[dt!=0])
  length(enhancers)
  cEnhan <- unlist(strsplit(enhancers,"-"))
  idx_chrR <- seq(1, length(cEnhan)-2, by = 3)
  idx_startR <- seq(2, length(cEnhan)-1, by = 3)
  idx_endR <- seq(3, length(cEnhan), by = 3)
  chrR <- cEnhan[idx_chrR]
  start.R <- as.numeric(cEnhan[idx_startR])
  end.R <- as.numeric(cEnhan[idx_endR])
  
  gr.enhancer.control[[cell]] <- GRanges(
    seqnames = Rle(chrR),
    ranges = IRanges(start.R, end = end.R))
  
  #unique enh
  dt <- control.uniqueEnh[,cell]
  enhancers <- names(dt[dt!=0])
  length(enhancers)
  cEnhan <- unlist(strsplit(enhancers,"-"))
  idx_chrR <- seq(1, length(cEnhan)-2, by = 3)
  idx_startR <- seq(2, length(cEnhan)-1, by = 3)
  idx_endR <- seq(3, length(cEnhan), by = 3)
  chrR <- cEnhan[idx_chrR]
  start.R <- as.numeric(cEnhan[idx_startR])
  end.R <- as.numeric(cEnhan[idx_endR])
  
  gr.enhancer.controlUnique[[cell]] <- GRanges(
    seqnames = Rle(chrR),
    ranges = IRanges(start.R, end = end.R))
}

#--cell types contain control
dt <- c("celltype","nb_ours","nb_anno","overlap_ours","overlap_anno","freq.overlap_ours","freq.overlap_anno")
for (cell in names(gr.enhancer.celltype)){
  if(cell %ni% "neuron"){
    annoteType <- cell_epiAnno[cell]
  }else{
    annoteType <- c("EX","INH")
  }
  gr.ours <- gr.enhancer.control[[cell]]
  gr.anno <- gr.epiAnno[gr.epiAnno$celltype %in% annoteType & gr.epiAnno$group %in% "Control"]
  
  o <- findOverlaps(gr.ours, gr.anno)
  o
  
  overlap.ours <- unique(gr.ours[queryHits(o)])
  overlap.anno <- unique(gr.anno[subjectHits(o)])
  
  freq.overlap.ours <- length(overlap.ours)/ length(gr.ours)
  freq.overlap.anno <- length(overlap.anno)/ length(gr.anno)
  
  
  saveLine <- c(cell, length(gr.ours), length(gr.anno), length(overlap.ours), length(overlap.anno), freq.overlap.ours, freq.overlap.anno)
  
  dt<-rbind(dt, saveLine)
  
}
colnames(dt) <- dt[1,]
rownames(dt) <- dt[,1]
dt<-dt[-1,-1]
dt 
overlap.celltype[["withSwarupNC_exactOverlap_controlEnh"]] <- dt

#--cell types contain control only
dt <- c("celltype","nb_ours","nb_anno","overlap_ours","overlap_anno","freq.overlap_ours","freq.overlap_anno")
for (cell in names(gr.enhancer.celltype)){
  if(cell %ni% "neuron"){
    annoteType <- cell_epiAnno[cell]
  }else{
    annoteType <- c("EX","INH")
  }
  gr.ours <- gr.enhancer.controlUnique[[cell]]
  gr.anno <- gr.epiAnno[gr.epiAnno$celltype %in% annoteType & gr.epiAnno$group %in% "Control"] # Swarup dataframe is unique (ignore 'joint')
  
  o <- findOverlaps(gr.ours, gr.anno)
  o
  
  overlap.ours <- unique(gr.ours[queryHits(o)])
  overlap.anno <- unique(gr.anno[subjectHits(o)])
  
  freq.overlap.ours <- length(overlap.ours)/ length(gr.ours)
  freq.overlap.anno <- length(overlap.anno)/ length(gr.anno)
  
  
  saveLine <- c(cell, length(gr.ours), length(gr.anno), length(overlap.ours), length(overlap.anno), freq.overlap.ours, freq.overlap.anno)
  
  dt<-rbind(dt, saveLine)
  
}
rownames(dt) <- dt[,1]
dt<-dt[,-1]
dt
overlap.celltype[["withSwarupNC_exactOverlap_controlEnhUnique"]] <- dt



#---------------------------------------------annotate 6) Bing Ren Science Paper
file_bing<-"/home/xhan/data_cellrangerOut/epiResource/Bing_etal_Science_cCREs.bed"
if(T){
  table_input <- read.csv(file_bing, sep = "\t", header = F)
  dim(table_input) #6728106 # ! more than 544k because one CRE can be shared in cell types!
  table_input[table_input$V4 %in% "cCRE_544727",] #cCRE_189997
  dim(unique(table_input[,c("V1","V2","V3","V4")])) #544729 4
  
  gr.enhancer.bing<- GRanges(
    seqnames = Rle(table_input$V1),
    ranges = IRanges(table_input$V2, end = table_input$V3))
  
  gr.enhancer.bing
  gr.enhancer.bing <- unique(gr.enhancer.bing)
  #check overlapping
  o2 <- findOverlaps(gr.enhancer.bing, gr.enhancer)
  nb_annoEnhan <- length(unique(gr.enhancer.bing))
  overlap_in_ourEn <- unique(gr.enhancer[subjectHits(o2)])
  nb_overlap <- length(overlap_in_ourEn)
  
  overlap_in_ourEn <- paste(seqnames(overlap_in_ourEn), overlap_in_ourEn@ranges, sep = "-")
  length(overlap_in_ourEn)

  name <- "BingScience"
  overlapList[[name]] <- overlap_in_ourEn 
  refEnh[[name]] <- nb_annoEnhan
}



##----------------------------------calculate for all resources: corresponding overlappedEnh distri in our cell types
dt.all_OverlapSets <- list()
for (resource in names(overlapList)){
  gr.overlap_in_ourEn <- overlapList[[resource]]
  tmp <- unlist(strsplit(gr.overlap_in_ourEn,split = "-"))
  idx_chr <- seq(1, length(tmp)-2, by = 3)
  idx_start <- seq(2, length(tmp)-1, by = 3)
  idx_end <- seq(3, length(tmp), by = 3)
  chr <- tmp[idx_chr]
  start<- as.numeric(tmp[idx_start])
  end<- as.numeric(tmp[idx_end])
  
  gr.overlap_in_ourEn <- GRanges(seqnames = Rle(chr), ranges = IRanges(start, end = end))

  
  for (cell in names(gr.enhancer.celltype)){
    gr.cell <- gr.enhancer.celltype[[cell]]
    m <- findMatches(gr.overlap_in_ourEn, gr.cell)
    cellExist <- rep(0, length(gr.overlap_in_ourEn))
    cellExist[queryHits(m)] <- 1
    
    values(gr.overlap_in_ourEn) <- cbind(values(gr.overlap_in_ourEn), DataFrame(cellExist))  
  }
  colnames(values(gr.overlap_in_ourEn)) <- names(gr.enhancer.celltype)
  gr.overlap_in_ourEn
  
  out <- NULL
  out2 <- NULL
  toCount <- data.frame(values(gr.overlap_in_ourEn)@listData)
  toCount <- cbind(toCount, sum=rowSums(toCount))
  unique <- toCount[toCount$sum==1,]
  for (cell in names(gr.enhancer.celltype)){
    dt <- values(gr.overlap_in_ourEn)[,cell]
    freq <- length(dt[dt!=0])/ length(gr.overlap_in_ourEn)
    out <- cbind(out, c(cell, length(dt[dt!=0]), freq))
    
    ## count those celltype unique
    dt2 <- unique[,cell]
    freq2 <- length(dt2[dt2!=0])/ length(gr.overlap_in_ourEn)
    out2 <- cbind(out2, c(paste(cell,"Unique", sep = ""), 
                        length(dt2[dt2!=0]), freq2))
    
  }
  out <- rbind(out,out2)
  rownames(out) <- rep(c("celltype", "nb_enh", "freq"),2)
  
  dt.all_OverlapSets[[resource]] <- out
}
dt.all_OverlapSets


##----------Overlap Number Summary
overlapEnh <- NULL
dt.summary <- NULL
for (i in names(overlapList)){
  overlapEnh <- c(overlapEnh, overlapList[[i]])
  dt.summary <- rbind(dt.summary, c(refEnh[[i]], length(overlapList[[i]])))
}
overlapEnh<-unique(overlapEnh)
length(overlapEnh)
length(overlapEnh) / sum(dt1[-1,]$Freq)

dt.summary <- rbind(dt.summary, c("/",length(overlapEnh)))
rownames(dt.summary) <- c(names(overlapList), "allResources")
dt.summary <- cbind(dt.summary, freq = sprintf("%.4f", as.numeric(dt.summary[,2])/sum(dt1[-1,]$Freq)))
colnames(dt.summary) <- c("nb_refEnh","nb_overlap", "freq_of_ourSet")
dt.summary


## save all this csv:
dt.summary
dt.all_OverlapSets
overlap.celltype

#write.csv(dt.summary, file.path(getwd(), projdir, "CRE/csv/validateCRE_summary.csv"))
#saveRDS(dt.all_OverlapSets, file.path(getwd(), projdir, "CRE/csv/validateCRE_overlapDistri.rds"))
#saveRDS(overlap.celltype, file.path(getwd(), projdir, "CRE/csv/validateCRE_exactOverlap.rds"))

###plot bar
setwd("/geschwindlabshares/RexachGroup/Xia_Data/atac_cellrangerOut/project0318")
projdir <- "projrmSubset_subPeak" 
dt.summary <- read.csv(file.path(getwd(), projdir, "CRE/csv/validateCRE_summary.csv"))

#223268 Enh
ordered_sources <- dt.summary[order(dt.summary$freq_of_ourSet),]$X

colors.condi <- c("darkred","lightgrey") 
names(colors.condi) <- c("Known","Novel")

dtp <- dt.summary[,c("X","freq_of_ourSet")]
colnames(dtp) <- c("source","Known")

dtp$Novel <- 1- dtp$Known
dtp<-reshape2::melt(dtp)
colnames(dtp) <- c("source","condi","freq")
dtp
dtp$source <- factor(dtp$source, levels=ordered_sources)
dtp$condi <- factor(dtp$condi, levels=c("Novel","Known"))

dtp$label <- dtp$freq * 100
dtp[dtp$condi %in% "Novel",]$label <- ""

p1 <- ggplot(data=dtp, aes(x=source, y = freq * 100, fill=condi)) + 
  geom_bar(stat="identity", position = "stack") +
  scale_fill_manual(values = colors.condi)  +
  labs(x="",y="",fill="",title="Percentage of Enhancers (%) #223268")+ 
  geom_text(aes(label=label), hjust = 0.5, vjust=-0.5, size = 2) +
  theme_set(theme_bw())+
  theme(plot.title = element_text(size=6),
        axis.title.x=element_text(size=6,face = "plain"),
        axis.title.y=element_text(size=6,vjust = 2,hjust = 0.5,face = "plain"),
        axis.text.x=element_text(size=6, face="plain", vjust = 1,hjust = 1, angle = 45),
        axis.text.y=element_text(size=6,face="plain"),
        legend.position = "right",
        #legend.position = c(1,1),
        #legend.title=element_blank(),
        legend.text=element_text(size=6),
        panel.grid =element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=1, colour = "black"))
pdf(file.path(file.path(getwd(), projdir, "bar_validateCRE.pdf")),
    width = 5.5, height = 2)
p1
dev.off()

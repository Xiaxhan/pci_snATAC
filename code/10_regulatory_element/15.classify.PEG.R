#1)PE types:    one/more to one/more relationshps for E-P P-G --> Five combined types
##  "df.PEG"
#2) PE existense:  majorgroups + subgroups
##  "dt.glCRE.MajorGroupExist" "dt.glCRE.groupExist"
#3) PE differential state :  P/E is DA? G is DEG?
## "PEGinfo"

library(ggplot2)
add<-theme_set(theme_bw())+
  theme(axis.title.x=element_text(size=13,face = "bold"),
        axis.title.y=element_text(size=13,vjust = 2,hjust = 0.5,face = "bold"),
        axis.text.x=element_text(size=13, face="bold", angle = 0, vjust=0.5, hjust = 0.5),
        axis.text.y=element_text(size=13,face="bold"),
        legend.position = "right",
        #legend.position = c(1,1),
        legend.title=element_blank(),
        legend.text=element_text(size=13),
        plot.title = element_text(size=15,face = "bold"),
        panel.grid =element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size=1, colour = "black"))
library(ArchR)
library(parallel)
addArchRThreads(threads = 10)
addArchRGenome("hg38")
library(ComplexHeatmap)
library(ggsignif)
library("RColorBrewer")

projdir <- "projrmSubset_subPeak"
proj <- loadArchRProject(projdir)

dt.glCRE.groupExist <- readRDS(file.path(getwd(), projdir, "CRE","subC.glCRE_groupExist.rds"))
dt.Enhancer.majorGroupExist <- readRDS(file.path(getwd(), projdir, "CRE","subC.enhancer_majorGroupExist.rds"))

dim(dt.Enhancer.majorGroupExist)

v.PE <- rownames(dt.glCRE.groupExist)
tmp <- unlist(strsplit(v.PE, ":"))
v.Promoter <-tmp[seq(1, length(tmp)-1, by=2)]
v.Eh <-tmp[seq(2, length(tmp), by=2)]
length(unique(v.Eh))

tmp <- unlist(strsplit(v.Promoter, "-"))
tmp.chr <-tmp[seq(1, length(tmp)-2, by=3)]
tmp.start <-tmp[seq(2, length(tmp)-1, by=3)]
tmp.end <-tmp[seq(3, length(tmp), by=3)]
gr.Promoter <- GRanges(seqnames = Rle(tmp.chr), ranges = IRanges(as.numeric(tmp.start), as.numeric(tmp.end)))
gr.Promoter

tmp <- unlist(strsplit(v.Eh, "-"))
tmp.chr <-tmp[seq(1, length(tmp)-2, by=3)]
tmp.start <-tmp[seq(2, length(tmp)-1, by=3)]
tmp.end <-tmp[seq(3, length(tmp), by=3)]
gr.Enh <- GRanges(seqnames = Rle(tmp.chr), ranges = IRanges(as.numeric(tmp.start), as.numeric(tmp.end)))
gr.Enh

# PEs are picked by P'overlappedG(if_exist)=NearestGene=expr_correlatedGene, 
# so the mathced peak's nearest gene is correct
m <- findMatches(gr.Promoter, proj@peakSet)
v.Genes <- proj@peakSet[subjectHits(m)]$nearestGene

df.PEG <- data.frame(PE = v.PE, P = v.Promoter, E = v.Eh, G = v.Genes)
#--------------------------------------------#1) PE types
#------1. check 1/more_To_1/more relationships
# P - G
tmp <- unique(data.frame(P=df.PEG$P, G=df.PEG$G))
#(as expected) 1P only match to 1G
count <- as.data.frame(table(tmp$P))
count[count$Freq>1,]
#1G could have >1 Ps
count <- as.data.frame(table(tmp$G))
count[count$Freq>1,]

# P - E
tmp <- unique(data.frame(P=df.PEG$P, E=df.PEG$E))
#1P regulated by >1E
count <- as.data.frame(table(tmp$P))
count[count$Freq>1,]
#1E regulate >1Ps
count <- as.data.frame(table(tmp$E))
count[count$Freq>1,]

#------2. add 'metainfo' for P/E/G types (for )
#1G could be regulated by >1 Ps
tmp <- unique(data.frame(P=df.PEG$P, G=df.PEG$G))
count <- as.data.frame(table(tmp$G))
GenebyP_Des <- count$Freq #Description
GenebyP_Des[GenebyP_Des>1] <- "m"
count <- cbind(count, GenebyP_Des=GenebyP_Des)
rownames(count) <- count$Var1

df.PEG<-cbind(df.PEG, GenebyP_NbP = count[df.PEG$G,]$Freq,  GenebyP_Des = count[df.PEG$G,]$GenebyP_Des)

# a P regulated by 1/m Es?
tmp <- unique(data.frame(P=df.PEG$P, E=df.PEG$E))
count <- as.data.frame(table(tmp$P))
PbyE_Des <- count$Freq #Description
PbyE_Des[PbyE_Des>1] <- "m"
count <- cbind(count, PbyE_Des=PbyE_Des)
rownames(count) <- count$Var1

df.PEG<-cbind(df.PEG, PbyE_NbE = count[df.PEG$P,]$Freq,  PbyE_Des = count[df.PEG$P,]$PbyE_Des)

# a E target(regulate) 1/m Ps?
count <- as.data.frame(table(tmp$E))
EtargetP_Des <- count$Freq #Description
EtargetP_Des[EtargetP_Des>1] <- "m"
count <- cbind(count, EtargetP_Des=EtargetP_Des)
rownames(count) <- count$Var1

df.PEG <- cbind(df.PEG, EtargetP_NbP = count[df.PEG$E,]$Freq,  EtargetP_Des = count[df.PEG$E,]$EtargetP_Des)

PEGtypes <- paste("GbyP_",df.PEG$GenebyP_Des, ":PbyE_",df.PEG$PbyE_Des, ":EtP_",df.PEG$EtargetP_Des, sep = "")
unique(PEGtypes)
df.PEG <- cbind(df.PEG, PEGtypes)

#--------------------------------------------#2) PE existence in major group
diag.flags <- c("Control" = 1000, "AD" = 100, "bvFTD" = 10, "PSP_S" = 1)
dt.diagnoseFlag<-NULL
for(i in 1:length(diag.flags)){
  symbol1 <- diag.flags[i]
  dt.diagnoseFlag <- rbind(dt.diagnoseFlag, data.frame(diagnose = names(symbol1), 
                                                       flag = symbol1))
}
for(i in 1:length(diag.flags)){
  symbol1 <- diag.flags[i]
  #print(paste(names(symbol1), symbol1, sep = ":"))
  s <- i+1
  if(i == length(diag.flags)){next}
  for (j in s:length(diag.flags)){
    symbol2 <- diag.flags[j]
    name <- paste(names(symbol1), names(symbol2), sep = "+")
    flag <- as.numeric(symbol1) + as.numeric(symbol2)
    dt.diagnoseFlag <- rbind(dt.diagnoseFlag, data.frame(diagnose = name, 
                                                         flag = flag))
  }
}
dt.diagnoseFlag <- rbind(dt.diagnoseFlag, data.frame(diagnose = c("Control+AD+bvFTD", "Control+bvFTD+PSP_S", 
                                                                  "AD+bvFTD+PSP_S","Control+AD+bvFTD+PSP_S"),
                                                     flag = c(1110, 1011, 111, 1111)))
dt.diagnoseFlag

head(dt.glCRE.groupExist)
dt.glCRE.MajorGroupExist <- NULL # label specified for diagnoses
j<-0
colvectorE <- NULL
for (major in c("ast", "mg", "neuron", "odc", "opc")){
  groups <- colnames(dt.glCRE.groupExist)[colnames(dt.glCRE.groupExist) %like% major & colnames(dt.glCRE.groupExist) %like% "-"]
  print(major)
  v_diagnoseExist <- rep(0, dim(dt.glCRE.groupExist)[1])
  dt <-  dt.glCRE.groupExist[,groups]
  for (condi in c("Control","AD","bvFTD","PSP_S")){
    dt.group <- dt[,colnames(dt) %like% condi]
    flag <- dt.diagnoseFlag[dt.diagnoseFlag$diagnose %in% condi,]$flag
    dt.group[dt.group==1] <- flag
    dt.group <- cbind(dt.group, sumflag=rowSums(dt.group))
    dt.group[,"sumflag"][dt.group[,"sumflag"] > 0] <- flag
    #unique(dt.group[,"sumflag"])
    v_diagnoseExist <- v_diagnoseExist + dt.group[,"sumflag"]
  }
  dt.glCRE.MajorGroupExist <- cbind(dt.glCRE.MajorGroupExist, v_diagnoseExist)
}
head(dt.glCRE.MajorGroupExist)
rownames(dt.glCRE.MajorGroupExist) <- rownames(dt.glCRE.groupExist)
colnames(dt.glCRE.MajorGroupExist) <- c("ast", "mg", "neuron", "odc", "opc")
dim(dt.glCRE.MajorGroupExist)
saveRDS(dt.glCRE.MajorGroupExist, file.path(getwd(), projdir, "CRE","subC.glCRE_MajorGroupExist.rds"))

dt.glCRE.MajorGroupExist <- readRDS(file.path(getwd(), projdir, "CRE","subC.glCRE_MajorGroupExist.rds"))

dt.glCRE.DxExist <- NULL
for (dx in c("AD","bvFTD","PSP_S","Control")){
  #dx <- "AD"
  groups <- colnames(dt.glCRE.groupExist)[colnames(dt.glCRE.groupExist) %like% dx & colnames(dt.glCRE.groupExist) %like% "-"]
  print(dx)
  dt.glCRE.groupExist[,groups][1:4,1:6]
  glCRE.dx <- dt.glCRE.groupExist[,groups]
  v.flag <- rowSums(glCRE.dx)
  v.flag[v.flag>0] <- 1
  dt.glCRE.DxExist <- cbind(dt.glCRE.DxExist, dx = v.flag)
}
colnames(dt.glCRE.DxExist) <-  c("AD","bvFTD","PSP_S","Control")

#count
tmp <- rownames(dt.glCRE.MajorGroupExist) %>% strsplit(.,":") %>% unlist(.)
V.Promoters <- tmp[seq(1, length(tmp), by = 2)]
V.Enhancers <- tmp[seq(2, length(tmp), by = 2)]
tmpadd <- dt.glCRE.MajorGroupExist
tmpadd[tmpadd!=0]<-1

V.PE <- rownames(dt.glCRE.MajorGroupExist)
#tmpadd <- cbind(tmpadd, sum=rowSums(tmpadd))
countEle <- "V.PE" #V.Promoters V.Enhancers

library(reshape2)
#tmp.df <- cbind(tmpadd, V.Promoters) %>% as.data.frame(.) #V.Promoters V.Enhancers
#tmp.df <- cbind(tmpadd, V.Enhancers) %>% as.data.frame(.)
tmp.df <- cbind(tmpadd, V.PE) %>% as.data.frame(.)
rownames(tmp.df) <- NULL
forcount <- melt(tmp.df, id.vars = countEle, variable.name="celltype",value.name="ifexist")
forcount <- forcount[forcount$ifexist!=0,] %>% unique(.)
table(forcount$celltype)
df.showtimes <- table(forcount[[countEle]]) %>% as.data.frame(.)
dim(df.showtimes[df.showtimes$Freq ==5,])[1] #common in 5 celltypes
unique <- df.showtimes[df.showtimes$Freq ==1,]$Var1
unique.df <- forcount[forcount[[countEle]] %in% unique,] 
table(unique.df$celltype)

#--------------------------------------------#3) PE differential state :  P/E is DA? G is DEG?
#----------snRNA DEG
deglist<-list()
deg_dir <- "/mypath/data_cellrangerOut/seurat_haplotype_dge"
for(region in c("preCG","insula")){
  for (celltype in c("astrocyte","excitatory","inhibitory","microglia","oligodendrocyte","opc")){
    filen <- paste(region, celltype, "haplotype.csv", sep="-")
    deglist[[region]][[celltype]] <- read.csv(file.path(deg_dir, filen), header = T)
  }
}

v.allGenes <- proj@geneAnnotation$genes$symbol
fdr_cutoff <- 0.1

v.deg<-NULL
dt.DEGanno_allGenes <- NULL
count.deg<-NULL
count.deg.filteredNAnames <- NULL
colnamestmp<-NULL
for(region in c("preCG","insula")){
  for (celltype in c("astrocyte","excitatory","inhibitory","microglia","oligodendrocyte","opc")){
    dt<-deglist[[region]][[celltype]]
    counts <- NULL
    counts2 <- NULL
    for(condi in c("clinical_dxAD","clinical_dxbvFTD","clinical_dxPSP.S","Tau_HH1.H2")){
      colnamestmp <- c(colnamestmp, paste(region,celltype, condi, sep = "-"))
      
      colname_fc <- paste(condi,".estimate",sep = "")
      colname_fdr <- paste(condi, ".p.value.adj", sep = "")
      dt.condi <- data.frame(gene = dt[["gene"]], log2FC=dt[[colname_fc]], fdr=dt[[colname_fdr]])

      vtmp <- rep("/", length(v.allGenes))
      names(vtmp) <- v.allGenes
      
      sigDEG <- dt.condi[dt.condi[["fdr"]] < fdr_cutoff,]
      v.deg <- c(v.deg, sigDEG$gene)
      sameGeneName <- intersect(sigDEG$gene, names(vtmp)) ### drop DEG name different from annotation
      sigDEG <- sigDEG[sigDEG$gene %in% sameGeneName,]
      
      vtmp[sigDEG$gene] <- sigDEG$log2FC
      dt.DEGanno_allGenes <- cbind(dt.DEGanno_allGenes, vtmp)
      counts <- cbind(counts, dim(dt.condi[dt.condi[["fdr"]] < fdr_cutoff,])[1]) # count the real DEG
      counts2 <- cbind(counts2, length(sigDEG$gene))# count DEGs with matched geneNames
    }
    count.deg <- rbind(count.deg, counts)
    count.deg.filteredNAnames <- rbind(count.deg.filteredNAnames, counts2)
  }
}
colnames(dt.DEGanno_allGenes) <- colnamestmp 
v.deg <-unique(v.deg)

rownames(count.deg) <- c(paste(c("ast","ex","inh","mg","odc","opc"), "-preCG", sep = ""),
                         paste(c("ast","ex","inh","mg","odc","opc"), "-insula", sep = ""))
colnames(count.deg) <- c("clinical_dxAD","clinical_dxbvFTD","clinical_dxPSP.S","Tau_HH1.H2")
count.deg
count.deg.filteredNAnames


#load sigDAR by region
sigDAR<-list()
sigDAR[["preCG"]] <- readRDS(file.path(getwd(), projdir,  "sigDAR.PreCG.rds"))
sigDAR[["midInsula"]] <- readRDS(file.path(getwd(), projdir,  "sigDAR.midInsula.rds"))

gr.allPeaks <- proj@peakSet #paste(seqnames(proj@peakSet), proj@peakSet@ranges)
gr.DARanno_allPeaks <- proj@peakSet
values(gr.DARanno_allPeaks) <- NULL
colnamesTmp <- NULL
for(region in c("preCG","midInsula")){
  for(group in names(sigDAR[[region]])){ #AD_ast.C1
    colnamesTmp <- c(colnamesTmp, paste(region, group, sep = "-"))
    gr <- sigDAR[[region]][[group]][[group]]
    m <- findMatches(gr, gr.allPeaks)
    m
    saveV <- rep("/", length(gr.DARanno_allPeaks))
    saveV[subjectHits(m)] <- gr$Log2FC 
    values(gr.DARanno_allPeaks) <- DataFrame(values(gr.DARanno_allPeaks), saveV)
  }
}
colnames(values(gr.DARanno_allPeaks)) <- colnamesTmp


## integrate to df.PEG
PEGinfo <- list()
PEGinfo[["dfPEG"]] <- df.PEG

m<-findMatches(gr.Promoter, gr.DARanno_allPeaks)
PEGinfo[["DA_P"]] <- values(gr.DARanno_allPeaks[subjectHits(m)])
rownames(PEGinfo[["DA_P"]]) <- paste(seqnames(gr.Promoter), gr.Promoter@ranges, sep="-")
m<-findMatches(gr.Enh, gr.DARanno_allPeaks)
PEGinfo[["DA_E"]] <- values(gr.DARanno_allPeaks[subjectHits(m)])
rownames(PEGinfo[["DA_E"]]) <- paste(seqnames(gr.Enh), gr.Enh@ranges, sep="-")

PEGinfo[["DA.gr"]] <- values(gr.DARanno_allPeaks)
rownames(PEGinfo[["DA.gr"]]) <- paste(seqnames(gr.DARanno_allPeaks), gr.DARanno_allPeaks@ranges, sep="-")

PEGinfo[["DEG"]] <- dt.DEGanno_allGenes[df.PEG$G,]

saveRDS(PEGinfo, file.path(getwd(), projdir, "CRE","subC.PEGtypes.rds"))

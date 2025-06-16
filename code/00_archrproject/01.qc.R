#Set QC cutoffs for individual samples

library(ArchR)
library(parallel)
addArchRThreads(threads = 10)
addArchRGenome("hg38")

outp <- "path_to_outdir"
setwd(outp)

# create new arrow files for samples
mysample <- read.csv("/home/xhan/atac/archrproject0318/list.qcsample", sep = "\t")
filename <- paste(mysample$sample, "fragments.tsv.gz", sep = ".")
frag_p<-"/home/xhan/data_cellrangerOut/fragments"
mysample <- data.frame(mysample, nfrags = 10^(mysample$log_nfrags))
mysample <- data.frame(mysample, path = paste(frag_p, filename, sep = "/"))
rownames(mysample) <- mysample$sample
mysample

ArrowFiles <- NULL
for (sp in mysample$sample){
  #sp <- "cat1_8"
  print(mysample[sp,]$path)
  ArrowFile <- createArrowFiles(
    inputFiles = mysample[sp,]$path,
    sampleNames = rownames(mysample[sp,]),
    minTSS =  mysample[sp,]$tssenrich,
    minFrags = 10^(mysample[sp,]$log_nfrags)
  )
  ArrowFiles <-c(ArrowFiles, ArrowFile)
}

# load ArchRproj and check QC plot
proj <- loadArchRProject("path_to_outdir/projALL")

df <- getCellColData(proj, select = c("log10(nFrags)", "TSSEnrichment", "Sample"))
ncell_df <-  dim(df)[1]
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 2, lty = "dashed") + 
  geom_vline(xintercept = 3, lty = "dashed") +
  ggtitle(paste("All 86 samples with", ncell_df, "cells", sep = " "))
plotPDF(p, name = "all.TSS-vs-Frags.pdf", ArchRProj = proj, addDOC = FALSE)

# plot for each condition
metainfo <- read.csv("~/atac/snATAC_metadata_summary2021_f.csv", header = T)
spinfo <- data.frame(sample = metainfo$ATAC_fastq_name, 
                     condition = metainfo$Clinical.Dx)
conditions <- unique(spinfo$condition)
for (condi in conditions){
  samples <- spinfo$sample[spinfo$condition %in% condi]
  samples
  
  dt <- df[df$Sample %in% samples,]
  n_cells <- dim(dt)[1]
  p <- ggPoint(
    x = dt[,1], 
    y = dt[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
  ) + geom_hline(yintercept = 2, lty = "dashed") + 
    geom_vline(xintercept = 3, lty = "dashed") +
    ggtitle(paste(length(samples), condi, "samples with", n_cells, "cells", sep = " "))
  fname <- paste(condi, "TSS-vs-Frags.pdf", sep = ".")
  plotPDF(p, name = fname, ArchRProj= proj, addDoc=F, width=5, height=5)
}

# plot for each sample
sample_set <- unique(df$Sample)
for (sample in sample_set){
  dt <- df[df$Sample %in% sample,]
  n_cells <- dim(dt)[1]
  p <- ggPoint(
    x = dt[,1], 
    y = dt[,2], 
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
    ylim = c(0, quantile(df[,2], probs = 0.99))
  ) + geom_hline(yintercept = 2, lty = "dashed") + 
    geom_vline(xintercept = 3, lty = "dashed") + 
    ggtitle(paste("Sample", sample, "with", n_cells, "cells", sep = " "))
  fname <- paste(sample, "TSS-vs-Frags.pdf", sep = ".")
  plotPDF(p, name = fname, ArchRProj= proj, addDoc=F, width=5, height=5)
}

 

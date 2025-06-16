#--Permutation test for "whether A disease peak set overrepresent/downrepresent methylated peaks (over - LOH; down - GOH)"

args <- commandArgs()
group <- args[6]
print(group)

library(dplyr)
library(ArchR)
library("stringr")
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsignif)
outp <- "sn_mC"

peakMC <- readRDS(file.path(outp, "peak_mClevel_hyperScore.rds"))

Groups <- unique(peakMC$group)
mC_refs <- unique(peakMC$mC_cell)
permutations <- 1:10000

length(Groups) * length(mC_refs) * 3


jid <- 0
dt.empiricalP <- NULL
for (refMC in mC_refs){
  dt.control <- peakMC[peakMC$mC_cell%in% refMC & peakMC$dx %in% "Control" & peakMC$group %in% group,]
  
  for (dx in c("AD","bvFTD","PSP_S")){
    dt.dx <- peakMC[peakMC$mC_cell%in% refMC & peakMC$dx %in% dx &peakMC$group %in% group,]
    if(!isEmpty(dt.dx[1])){
      jid <- jid + 1
      # make dataframe for the original two groups - counts/ freq/ log2(freq1/freq2)
      counts <- rbind(table(dt.dx$label_hyper) %>% as.data.frame(.),
                      table(dt.control$label_hyper) %>% as.data.frame(.))
      counts$Var1 <- c("treat_WOmc","treat_Wmc","ct_WOmc", "ct_Wmc")
      counts <- rbind(counts, c("ct_sum", dim(dt.control)[1]), c("treat_sum", dim(dt.dx)[1]))
      colnames(counts) <- c("type","original")
      rownames(counts) <- counts$type
      treat_freq <- as.numeric(counts["treat_Wmc",]$original) / as.numeric(counts["treat_sum",]$original)
      ct_freq <- as.numeric(counts["ct_Wmc",]$original) / as.numeric(counts["ct_sum",]$original)
      #fc <- treat_freq/ct_freq #
      freq_diff <- treat_freq - ct_freq
      counts <- rbind(counts, 
                      c("treat_freq", treat_freq),
                      c("ct_freq", ct_freq),
                      c("freq_diff", freq_diff))
      rownames(counts) <- counts$type
      counts
      
      
      dt.pool <- rbind(dt.control, dt.dx)
      dt.All <- counts$original %>% as.numeric() %>% as.data.frame() %>% t(.)
      rownames(dt.All) <- "original"
      colnames(dt.All) <- rownames(counts)
      dt.All <- as.data.frame(dt.All)
      dt.All
      
      # permutation to get mC-peak counts -> only nb
      TreatSize <- dt.All["original","treat_sum"]
      sampleSize <- dim(dt.pool)[1]
      dt.permutation <- NULL
      for (i in permutations){
        set.seed(i)
        print(paste("jid:",jid, "---permutation:",i,"---",group, ".",refMC, ".", dx, sep = ""))
        # permute to create a new "treat" and the remaining are "control"
        perm_treat <- sample(sampleSize, TreatSize) 
        perm_ct <- setdiff(1:sampleSize, perm_treat) 
        count1 <- dt.pool$label_hyper[perm_treat] %>% table(.) %>% as.data.frame(.)
        count2 <- dt.pool$label_hyper[perm_ct] %>% table(.) %>% as.data.frame(.)
        treat_Wmc <- count1[count1$. %in% 1,]$Freq
        ct_Wmc <- count2[count2$. %in% 1,]$Freq
        dt.permutation <- rbind(dt.permutation, c(treat_Wmc, ct_Wmc))
      }
      
      rownames(dt.permutation) <- paste("permute",permutations, sep = "")
      colnames(dt.permutation) <- c("treat_Wmc", "ct_Wmc")
      dt.permutation
      dt.permutation <- as.data.frame(dt.permutation)
      
      # calculate the freq and "measurement"
      treat_freq <- dt.permutation$treat_Wmc / dt.All["original","treat_sum"]
      ct_freq <- dt.permutation$ct_Wmc / dt.All["original","ct_sum"]
      #fc <- treat_freq/ct_freq 
      #quantile(dt.permutation$freq_diff)
      freq_diff <- treat_freq - ct_freq  
      quantile(freq_diff)
      
      dt.permutation <- cbind(dt.permutation, treat_freq, ct_freq, freq_diff)
      dt.permutation
      
      # calculate the empirical P for each group
      right.p <- length(freq_diff[freq_diff >= dt.All$freq_diff]) / length(permutations)
      left.p <- length(freq_diff[freq_diff <= dt.All$freq_diff]) / length(permutations)
      
      dt.empiricalP <- rbind(dt.empiricalP, 
                             data.frame(group, refMC, dx, left.p, right.p, 
                                        freqdiff_dx2ct=dt.All$freq_diff, 
                                        percent95=quantile(freq_diff, 0.95),
                                        percent5=quantile(freq_diff, 0.05) ) )
    }
  }
  
}

fname <- paste("permutation10k", group, "csv", sep = "_")
write.csv(dt.empiricalP, row.names = F, file.path(outp, fname))
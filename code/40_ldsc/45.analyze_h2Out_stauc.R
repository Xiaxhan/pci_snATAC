library(stringr)
library(data.table)
library(ggplot2)
library(ggpubr)
outp <- "/geschwindlabshares/RexachGroup/Xia_Data/ldsc/data/plot"
nbSNP <- read.csv("/geschwindlabshares/RexachGroup/Xia_Data/ldsc/MS_more_dx/SNPnb_from_LDSC.txt", 
                  header = T, sep = "\t")
GWASs <- c("ADHD_Demontis_ng2023","ASD_Grove_ng2019", "MDD_Howard_nn2019", "SCZ_Trubetskoy_n2022",
           "MS_Andlauer_SciAdv2016", "ALS_Rheenen_ng2021", "IBD_Glanville_2021",
           "AD_Bellenguez_ng2022", "AD_Jansen", "AD_Kuncle", "AD_Lambert", "AD_MorenoGrau",  "PSP_Hoglinger",
           "FTD_META")
           #"FTD_bvFTD", "FTD_MMD", "FTD_SD", "FTD_PNFA", )
mydir <- "/geschwindlabshares/RexachGroup/Xia_Data/ldsc/MS_more_dx"
#invariant peaks:
list1 <- c("list.majorC","list.50subC") #list.50subC -> re run for invariant P
inp1 <- "/geschwindlabshares/RexachGroup/Xia_Data/ldsc/data/outInvariantP_h2"
#dynamic peaks:
list2 <- c("list.majorC","list.50subC","list.majorCupdown","list.majorCdx")
inp2 <- "/geschwindlabshares/RexachGroup/Xia_Data/ldsc/data/outMarkerP_h2"

h2_result <- list()
#Invariant
for (file in list1){
  type <- gsub("list.","",file)
  groups <- read.csv(file.path(mydir, file), header = F)$V1
  for (group in groups){
    for (gwas in GWASs){
      ldsc_file <- paste(group, gwas, "results", sep = ".")
      # subC has noresults for "AD_MorenoGrau.results"
      if(file.exists(file.path(inp1, ldsc_file))){
        dt <-  read.csv(file.path(inp1, ldsc_file), sep = "\t")
        group_gwas <- paste(type, "InvariantP",group, gwas, sep = "|")
        #print(group_gwas)
        h2_result[[group_gwas]] <- dt
      }else{
        print(paste("No file:",ldsc_file))
      }

    }
  }
}

#Dynamic
for (file in list2){
  type <- gsub("list.","",file)
  groups <- read.csv(file.path(mydir, file), header = F)$V1
  for (group in groups){
    for (gwas in GWASs){
      ldsc_file <- paste(group, gwas, "results", sep = ".")
      # subC has noresults for "AD_MorenoGrau.results"
      if(file.exists(file.path(inp2, ldsc_file))){
        dt <-  read.csv(file.path(inp2, ldsc_file), sep = "\t")
        group_gwas <- paste(type, "DynamicP",group, gwas, sep = "|")
        #print(group_gwas)
        h2_result[[group_gwas]] <- dt
      }else{
        print(paste("No file:",ldsc_file))
      }
      
    }
  }
}

#----
allC <- NULL ## L2_0! + 53 baseline-categories
df.h2 <- data.frame()
for (group_gwas in names(h2_result)){
  allC <- c(allC, h2_result[[group_gwas]]$Category)
  dt <- h2_result[[group_gwas]][1,]
  tmp <- unlist(strsplit(group_gwas, "\\|"))
  group_type <- tmp[1]
  PeakType <- tmp[2]
  group  <- tmp[3]
  gwas  <- tmp[4]
  # calculate *tau_c (the standardized tau_c)
  # enrichment=tau_c(C) / tau_c(overall)
  # the"C" - SNP_in_C: uses the total SNPs that LDSC uses, and x the props_of_SNPs
  enrich_score <- dt[,"Enrichment"]
  tau_c <-  dt[,"Coefficient"]
  tau_c_overall <- tau_c/ enrich_score
  Prop_SNPs <- dt[,"Prop._SNPs"]
  nb_SNPs <- nbSNP[nbSNP$GWAS %in% gwas,]$nSNPs_LDSC_read
  SNP_in_C <- c(rep(1, ceiling(Prop_SNPs*nb_SNPs)),
                rep(0, nb_SNPs - ceiling(Prop_SNPs*nb_SNPs)) )
  standardized_tau_c <-  (tau_c * sd(SNP_in_C)) / tau_c_overall
  
  df.h2 <- rbind(df.h2, data.frame(
    group_type, PeakType, group,
    GWAS = gwas,
    Enrichment = dt$Enrichment, 
    Enrichment_p = dt$Enrichment_p,
    Coefficient = dt$Coefficient,
    Coefficient_std_error = dt$Coefficient_std_error,
    Coefficient_zscore=dt$"Coefficient_z.score",
    standardized_tau_c) )
}
tail(df.h2)

GWASs <- unique(df.h2$GWAS)
select_GWASs <- c(GWASs[GWASs %like% "AD_"], "FTD_META", "PSP_Hoglinger")

# Split by Groups to conduct FDR for s_tau_c
df_majorC <- df.h2[df.h2$group_type %in% "majorC" & df.h2$GWAS %in% select_GWASs,]

# Calculate Pval for s_tau_c
if(T){
  df_used <- df_majorC
  
  v_stauc <- df_used$standardized_tau_c
  sqrt(sum((v_stauc - mean(v_stauc)) ^ 2/(length(v_stauc) - 1))) /sqrt(length(v_stauc))
  plotrix::std.error(v_stauc)
  se_stauc <- plotrix::std.error(v_stauc)
  stauc_variable <- v_stauc/ se_stauc
  stauc_p <- data.frame(v_stauc, stauc_variable, Pval=pnorm(-abs(stauc_variable)))
  #*tau_c fdr
  stau_c_fdr <- p.adjust(stauc_p$Pval, method = "BH", n = length(stauc_p$Pval))

  #  # *< 0.05; ** <0.005; *** < 0.001
  FDRs <- stau_c_fdr
  FDRs[FDRs<0.001] <- "***"
  FDRs[!FDRs %like% "\\*"][FDRs[!FDRs %like% "\\*"] < 0.005] <- "**"
  FDRs[!FDRs %like% "\\*"][FDRs[!FDRs %like% "\\*"] < 0.05] <- "*"
  FDRs[!FDRs %like% "\\*"][FDRs[!FDRs %like% "\\*"] > 0.05] <- ""
  
  #enrich-label
  pvalue <- df_used$Enrichment_p
  pvalue[pvalue<0.001] <- "***"
  pvalue[!pvalue %like% "\\*"][pvalue[!pvalue %like% "\\*"] < 0.005] <- "**"
  pvalue[!pvalue %like% "\\*"][pvalue[!pvalue %like% "\\*"] < 0.05] <- "*"
  pvalue[!pvalue %like% "\\*"][pvalue[!pvalue %like% "\\*"] > 0.05] <- ""
}
df_majorC <- cbind(df_majorC, 
                   Enrichment_p_flag = pvalue,
                   standardized_tau_c_Pval=stauc_p$Pval,
                   standardized_tau_c_fdr=stau_c_fdr,
                   stauc_fdr_flag = FDRs)

#-plot
#*tau_c/se(*tauc) distribution
# Normality test
shapiro.test(stauc_variable)
shapiro.test(rnorm(50, mean = 0, sd = 1))  
shapiro.test(rnorm(50, mean = 4, sd = 11)) 
df <- data.frame(values = stauc_variable)
p1 <- ggdensity(stauc_variable, 
                main = "Density plot",
                xlab = "*tau_c/sd(*tau_c)")
p2 <- ggplot(df, aes(x = values)) +
  geom_histogram(aes(y=..density..),
                 binwidth = 0.1, 
                 fill = "lightblue", 
                 color = "black",
                 size = 0.1,
                 alpha = 0.7) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), color = "red", size = 1) +
  labs(title = "majorC: Density histogram of *tau_c/se(*tau_c)", x = "*tau_c/sd(*tau_c)") +
  theme_minimal()
outfile <- "histogram.tauc_variable.majorC.pdf"
pdf(file.path(outp, outfile), width = 5, height = 3)
p1
p2
dev.off()

add<-theme_set(theme_bw())+
  theme(plot.title = element_text(size = 10, face = "bold"),
        axis.title.x=element_text(size=10,face = "plain"),
        axis.title.y=element_text(size=10,vjust = 2, hjust = 0.5,face = "plain"),
        axis.text.x=element_text(size=10, face="plain", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y=element_text(size=10,face="plain"),
        legend.position = "right",
        #legend.position = c(1,1),
        #legend.title=element_blank(),
        legend.text=element_text(size=10),
        #panel.grid =element_blank(),
        #panel.border = element_blank(),
        panel.border = element_rect(color = "gray"),
        axis.line = element_line(size=1, colour = "black"))
df_majorC$GWAS <- factor(df_majorC$GWAS, levels = sort(unique(df_majorC$GWAS)) ) 
df_majorC$group <- factor(df_majorC$group, levels = unique(df_majorC$group))
df_majorC$PeakType <- factor(df_majorC$PeakType, levels = unique(df_majorC$PeakType))

outfile <- "heamtaps.LDSC.majorC_dynamic_invariant.pdf"
titlename1 <- "LDSC GWAS enrichent in majorC dyanmic/invariant peaks, Pval < 0.05"
p1 <- ggplot(df_majorC, aes(x = group, y = GWAS, fill = Enrichment)) +
  geom_tile(colour = "black") + 
  geom_text(aes(label = Enrichment_p_flag)) +
  facet_grid(~PeakType) + 
  scale_fill_gradient2(midpoint = 0, low = "navy",  mid = "white", high = "firebrick3") +
  #geom_vline(xintercept = c(3.5, 8.5, 13.5, 18.5, 23.5, 28.5)) + 
  labs(fill = "Enrichment", title = titlename1) + add 

titlename2 <- "LDSC standardized tau_c in dyanmic/invariant peaks, FDR of *tau_c < 0.05"
p2 <- ggplot(df_majorC, aes(x = group, y = GWAS, fill = standardized_tau_c)) +
  geom_tile(colour = "black") + 
  geom_text(aes(label = stauc_fdr_flag)) +
  facet_grid(~PeakType) + 
  scale_fill_gradient2(midpoint = 0, low = "navy",  mid = "white", high = "firebrick3") +
  labs(fill = "*tau_c", title = titlename2) + add 

pdf(file.path(outp, outfile), width = 8, height = 3.5)
p1
p2
dev.off()

outfile2 <- "bar.LDSC.majorC_dynamic_invariant.pdf"
titlename2 <- "LDSC standardized tau_c in dyanmic/invariant peaks, τ* FDR: *< 0.05; ** <0.005; *** < 0.001"
df_majorC$standardized_tau_c_rmNegative <- df_majorC$standardized_tau_c
df_majorC[df_majorC$standardized_tau_c_rmNegative<0,]$standardized_tau_c_rmNegative <- 0
df_majorC[df_majorC$standardized_tau_c_rmNegative==0,]$stauc_fdr_flag <- ""

head(df_majorC)
mycolors <- c("#4857a7","grey")
names(mycolors) <- c("DynamicP","InvariantP")
df_majorC$PeakType <- factor(df_majorC$PeakType ,levels= rev(c("DynamicP","InvariantP")) )

p_bar <- ggplot(df_majorC, aes(y = standardized_tau_c_rmNegative, 
                            x = GWAS, 
                            fill = PeakType)) +
  geom_bar(stat="identity", position = "dodge")+ 
  geom_text(aes(label = stauc_fdr_flag), size = 3.3, position = position_dodge(width = .9)) +
  facet_grid(~group) + 
  scale_fill_manual(values = mycolors) +
  labs(y = "h^2 effect size (τ*)", title = titlename2) + add
pdf(file.path(outp, outfile2), width = 16, height = 3)
p_bar
dev.off()
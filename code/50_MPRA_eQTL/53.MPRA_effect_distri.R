# Distribution of MPRA snps
library(ggplot2)
library(dplyr)
library(data.table)
outp <- "QTL_fsVars/out"
p1 <- "resources"
f_snp <- "MPRA_list_of_variants.tsv"
f_effect <- "HMC3_MPRA_rep1rep2.tsv"
mpra_pos <- read.csv(file.path(p1,"MPRA_w_pos_GRCh38.csv")) # from 'MPRA_cal_fishertest"
dyn_inv_peaks_w_mpra <- read.csv(file.path(outp, "dynamicP_invariantP_w_MPRA.csv"), row.names=1)

df_snps <- read.csv(file.path(p1, f_snp), row.names = 1) #11504
df_effect <- read.csv(file.path(p1, f_effect), sep = "\t")
df_effect$lib <- gsub(".*==","",df_effect$Names)
rownames(df_effect) <- df_effect$lib
rownames(mpra_pos) <- mpra_pos$refsnp_id

library_undefine <- setdiff(df_snps$library,df_effect$lib)
library_shared <- intersect(df_snps$library,df_effect$lib)


df_snps <- df_snps[df_snps$library %in% library_shared, ] #11165
df_snps <- cbind(df_snps, 
                 df_effect[df_snps$library, ],
                 mpra_pos[df_snps$variant_id, c("chr_name","chrom_start")])
dim(unique(df_snps[,c("library","lib")]))

all(df_snps$library == df_snps$lib) 

colnames(df_snps)[1:3] <- c("chr_hg19","start_hg19","end_hg19")
df_snps$lib <- NULL
colnames(df_snps)[8] <- "Name_repLibrary"
colnames(df_snps)[12:13] <- c("chr_hg38","start_hg38") 

#Plot histograms
add <- theme_set(theme_bw())+
  theme( plot.title = element_text(size=12, face = "bold"),
         axis.title.x=element_text(size=12,face = "plain"),
         title =element_text(size=12, face='bold'),
         axis.title.y=element_text(size=12,vjust = 2,hjust = 0.5,face = "plain"),
         axis.text.x=element_text(size=12, face="plain", angle = 0, hjust = 0.5),
         #axis.text.x=element_blank(),
         axis.text.y=element_text(size=12,face="plain"),
         legend.position = "right",
         legend.title=element_blank(),
         legend.text=element_text(size=9),
         panel.grid = element_blank(),
         #panel.grid.minor.y = element_line(linewidth=1, colour = "grey"),
         #panel.border = element_rect(linetype = "dashed", fill = NA),
         panel.border = element_blank(),
         axis.ticks = element_line(colour = "black", linewidth= 0.1),
         axis.line = element_line(linewidth= 0.1, colour = "black", linetype="solid"))
p1 <- ggplot(df_snps, aes(x = pval)) +
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept=0.05, color = "red", linetype = "dashed") +
  labs(y=paste("Nb of variants(total:", dim(df_snps)[1])) + add
p2 <- ggplot(df_snps, aes(x = fdr)) +
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept=0.05, color = "red", linetype = "dashed") +
  labs(y="Nb of variants") + add
p3 <- ggplot(df_snps, aes(x = logFC)) +
  geom_histogram(binwidth = 0.01) + 
  geom_vline(xintercept=1, color = "red", linetype = "dashed") +
  geom_vline(xintercept=-1, color = "red", linetype = "dashed") +
  labs(y="Nb of variants") + add

pdf(file.path(outp, "histogram.MPRA.pdf"), width = 6, height = 3)
p1
p2
p3
dev.off()

#
dup_snps <- table(dyn_inv_peaks_w_mpra$snp)[table(dyn_inv_peaks_w_mpra$snp) > 1]
dyn_inv_peaks_w_mpra[dyn_inv_peaks_w_mpra$snp %in% names(dup_snps),]
dyn_inv_peaks_w_mpra <- dyn_inv_peaks_w_mpra[!(dyn_inv_peaks_w_mpra$snp %in% names(dup_snps) &
					       dyn_inv_peaks_w_mpra$label %in% "Invariant"),]

snps_in_dyn <- dyn_inv_peaks_w_mpra[dyn_inv_peaks_w_mpra$label %in% "Dynamic",]$snp %>%
	strsplit(.,split=";") %>% unlist() %>% sort() %>% unique()
snps_in_dyn <- intersect(snps_in_dyn, df_snps$variant_id)
snps_in_inv <- dyn_inv_peaks_w_mpra[dyn_inv_peaks_w_mpra$label %in% "Invariant",]$snp %>%
        strsplit(.,split=";") %>% unlist() %>% sort() %>% unique()
snps_in_inv <- intersect(snps_in_inv, df_snps$variant_id)
print(length(snps_in_dyn))
print(length(snps_in_inv))

plot_histogram <- function(label){
mpra_in_peak <- dyn_inv_peaks_w_mpra[dyn_inv_peaks_w_mpra$label %in% label,]$snp
mpra_in_peak <- unlist(strsplit(mpra_in_peak, split = ";"))
mpra_in_peak <- unique(mpra_in_peak)
mpra_in_peak <- intersect(mpra_in_peak, df_snps$variant_id)
df_snps_inP <- df_snps[df_snps$variant_id %in% mpra_in_peak,]
ylab <- paste("#variants in", label, "peaks", length(mpra_in_peak))
p1 <- ggplot(df_snps_inP, aes(x = pval)) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept=0.05, color = "red", linetype = "dashed") +
  labs(y=ylab) + add
p2 <- ggplot(df_snps_inP, aes(x = fdr)) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept=0.05, color = "red", linetype = "dashed") +
  labs(y=ylab) + add
p3 <- ggplot(df_snps_inP, aes(x = logFC)) +
  geom_histogram(binwidth = 0.01) +
  geom_vline(xintercept=1, color = "red", linetype = "dashed") +
  geom_vline(xintercept=-1, color = "red", linetype = "dashed") +
  labs(y=ylab) + add
fout <- paste("histogram.MPRA.in", label, "peaks", "pdf", sep = ".")
pdf(file.path(outp, fout), width = 6, height = 3)
print(p1)
print(p2)
print(p3)
dev.off()
}

plot_histogram("Dynamic")
plot_histogram("Invariant")

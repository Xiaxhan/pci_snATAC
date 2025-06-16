library(dplyr)
library(data.table)
library(GenomicRanges)
library(ggplot2)
outp <- "QTL_fsVars/out"
majorCs <- c("ast","mg","EX","IN","odc","opc")
Dxs <- c("Control","AD","bvFTD","PSP_S")
p1 <- "resources"
colors.peakType <- ArchR::ArchRPalettes$summerNight[1:5]
names(colors.peakType) <- c("Distal","Enhancer","Promoter","Exonic","Intronic")

mpra_in_celltypePeaks <- read.csv(file.path(outp, "mpra_in_celltypePeaks.csv"), row.names=1)
mpra_in_celltypePeaks <- mpra_in_celltypePeaks[!is.na(mpra_in_celltypePeaks$chr_hg38),]

fsVar <- mpra_in_celltypePeaks[mpra_in_celltypePeaks$pval < 0.05,] 

# Cell type specific peaks contain MPRA; split by peakType
dtp <- table(dtt$ct, dtt$peakType) %>% data.frame()
colnames(dtp) <- c("group","peakType","nb")
dtp$peakType <- factor(as.character(dtp$peakType), levels = names(colors.peakType))

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

p <- ggplot(dtp, aes(x = peakType, y=nb, fill=peakType)) +
  geom_bar(stat="identity", position="dodge", size=0.1) + 
  coord_flip() +
  facet_grid(~group) +
  labs(x="",y="",fill="", title="#Peaks contain MPRA frVar")+
  scale_fill_manual(values=colors.peakType) + add

pdf(file.path(outp, "bar.peaks_w_mpra_fsVar.pdf"), width = 10, height = 3)
p
dev.off()



# Cell type specific dynamic/invariant peaks with MPRA fsVar
dt2 <- fsVar[,c("loci","ct","label")] %>% unique()
table(dt2$ct, dt2$label)

#true enhancer not CRE
dt2_enh <- fsVar[fsVar$peakType %in% "Enhancer", c("loci","ct","label")] %>% unique()
table(dt2_enh$ct, dt2_enh$label)

# plot- peaks
get_dtp <- function(dt){
  dtp <- table(dt$ct, dt$label) %>% data.frame()
  colnames(dtp) <- c("group","peakType","nb")
  dtp$peakType <- factor(as.character(dtp$peakType), levels = c("Dynamic","Invariant"))
  return(dtp)
}
# plot-enh
colors.dyIn <- ArchR::ArchRPalettes$stallion2[1:2]
names(colors.dyIn) <- c("Dynamic","Invariant")
dtp_p <- get_dtp(dt2)
dtp_e <-  get_dtp(dt2_enh)
p2p <- ggplot(dtp_p, aes(x = group, y=nb, fill=peakType)) +
  geom_bar(stat="identity", position="dodge", size=0.1) +
  labs(x="",y="",fill="", title="#Peaks contain MPRA frVar")+
  scale_fill_manual(values=colors.dyIn) + add
p2e <- ggplot(dtp_e, aes(x = group, y=nb, fill=peakType)) +
  geom_bar(stat="identity", position="dodge", size=0.1) +
  labs(x="",y="",fill="", title="#Enhancer contain MPRA frVar")+
  scale_fill_manual(values=colors.dyIn) + add
f2p <- ggplot(dtp_p, aes(x = group, y=nb, fill=peakType)) +
  geom_bar(stat="identity", position="fill", size=0.1) +
  labs(x="",y="",fill="", title="%Peaks contain MPRA frVar")+
  scale_fill_manual(values=colors.dyIn) + add
f2e <- ggplot(dtp_e, aes(x = group, y=nb, fill=peakType)) +
  geom_bar(stat="identity", position="fill", size=0.1) +
  labs(x="",y="",fill="", title="%Enhancer contain MPRA frVar")+
  scale_fill_manual(values=colors.dyIn) + add

pdf(file.path(outp, "bar.Peak-ct_Enh_w_mpra_fsvar.pdf"), width = 5, height = 2)
p2p
p2e
f2p
f2e
dev.off()

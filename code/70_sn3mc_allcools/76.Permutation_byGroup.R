library(dplyr)
library(ArchR)
library("stringr")
library(ggplot2)
library(ggpubr)
library(reshape2)
library(ggsignif)
outp <- "sn_mC"
colors.dx <- ArchR::ArchRPalettes$circus[1:4]
names(colors.dx) <- c("Control","AD","bvFTD","PSP_S")

#---load 50 subCs
projmeta <- readRDS("projATAC_50subC_meta_584904.rds")
my_subcs <- unique(projmeta$subClusters) %>% sort()

#----------------------load in loop Run output
Groups <- read.csv(file.path(outp, "celltypes.csv"), header = F)$V1
# delete empty files
drop_files <- read.csv(file.path(outp, "drop_files.txt"), header = F)$V1

Groups <- c(Groups[Groups %like% "ast.C"],Groups[Groups %like% "mg.C"])
Groups <- Groups[Groups %in% my_subcs]

####  use subC-unique peak 
dt.empiricalP <- NULL
for (group in Groups){
  for (dx in c("AD","bvFTD","PSP_S")){
    # use subC-unique peak or not
    fname <- paste("permutation10k_subCunique", dx, group, "csv", sep = ".")
    
    if(fname %in% drop_files){
      next
    }else{
      
      if(file.exists(file.path(outp, fname))){
        dt <- read.csv(file.path(outp, fname), header = T)
        dt.empiricalP <- rbind(dt.empiricalP, dt)
      }else{
        print(fname)
      }
      
      
    }
  }
}

dt.empiricalP <- unique(dt.empiricalP)
dim(unique(dt.empiricalP)) 
table(dt.empiricalP$dx, dt.empiricalP$group)  

# 2) heatmap to show the siginifcant LOH/ROH
permute.left <- dt.empiricalP$left.p
permute.right <- dt.empiricalP$right.p

permute.left[permute.left <=  0.001] <- "***"
permute.left[permute.left > 0.001 & permute.left <= 0.01] <- "**"
permute.left[permute.left > 0.01 & permute.left <= 0.05] <- "*"
permute.left[!permute.left %like% "\\*"] <- NA

permute.right[permute.right <=  0.001] <- "***"
permute.right[permute.right > 0.001 & permute.right <= 0.01] <- "**"
permute.right[permute.right > 0.01 & permute.right <= 0.05] <- "*"
permute.right[!permute.right %like% "\\*"] <- NA

v.p <- 1:length(permute.left)

dt$max_pval <- pmin(dt$left.p, dt$right.p)
dt$fdr <- p.adjust(dt$max_pval, method = "fdr")

dt[dt$group %in% "mg.C4" & dt$refMC %like% "Micro", ]
dt[dt$group %in% "ast.C1" & dt$refMC %like% "Astro", ]

dt$fdr_label <- ifelse(dt$fdr<0.01, "*","")

write.csv(dt, row.names = F, 
          file.path(outp, "uniqueMC.ASC_MG.DXpeak_mCfreq.permutation.csv"))

# plot 
dt$showP <- ifelse(dt$P %in% "***", "*", "")
dt$label <- ifelse(dt$freqdiff_dx2ct > 0, "Loss", "Gain")
dt$abs_freqdiff_dx2ct <- abs(dt$freqdiff_dx2ct) * 100

dt_asc <- dt[dt$group %like% "ast" & dt$refMC %in% "Astro", ]
dt_mg <- dt[dt$group %like% "mg" & dt$refMC %in% "NonN_Micro-Endo_TYROBP", ]
dt_asc$group <- factor(dt_asc$group , levels=str_sort( unique(dt_asc$group) , numeric = TRUE) )
dt_mg$group <- factor(dt_mg$group , levels=str_sort( unique(dt_mg$group) , numeric = TRUE) )


add <- theme_set(theme_bw())+
  theme( plot.title = element_text(size=9, face = "plain"),
         axis.title.x=element_text(size=9,face = "plain"),
         title =element_text(size=9, face='bold'),
         axis.title.y=element_text(size=9,vjust = 2,hjust = 0.5,face = "plain"),
         axis.text.x=element_text(size=9, face="plain", angle = 90, hjust = 0.5),
         #axis.text.x=element_blank(),
         axis.text.y=element_text(size=9,face="plain"),
         legend.position = "right",
         
         legend.title=element_blank(),
         legend.text=element_text(size=9),
         panel.grid = element_blank(),
         #panel.grid.minor.y = element_line(linewidth=1, colour = "grey"),
         #panel.border = element_rect(linetype = "dashed", fill = NA),
         panel.border = element_blank(),
         axis.ticks = element_line(colour = "black", linewidth= 0.1),
         axis.line = element_line(linewidth= 0.1, colour = "black", linetype="solid"))


colors.LOH <- c("#e51e1b","#387fb8")
names(colors.LOH) <- c("Loss","Gain")

max(dt_asc$abs_freqdiff_dx2ct )
max(dt_mg$abs_freqdiff_dx2ct )

p_asc <- list()
for (dx in c("AD","bvFTD","PSP_S")){
  dtss <- dt_asc[dt_asc$dx %in% dx,]
  p_asc[[dx]] <- ggplot(dtss, aes(x = group, y = abs_freqdiff_dx2ct, fill = label)) +
    geom_bar(stat="identity", position = "dodge", color = NA) +
    ylim(0,1) +
    labs(x = "", y = "Abs. Diff. in mC Peak %", 
         title = paste(dx, ": Astrocyte-specific heterochromatin gain or loss 
       (empirical P <= 0.001 by Permutation 10k )", sep ="" )) +
    geom_text(# position = position_dodge(width = 0.7),  # Adjust text position relative to bars
              size = 5, hjust=0.5,
              #aes(label = fdr_label) ) +
              #aes(label = showP) ) +
              aes(label = P) ) +
    scale_fill_manual(values = colors.LOH) + add
  
}

p_mg <- list()
for (dx in c("AD","bvFTD","PSP_S")){
  dtss <- dt_mg[dt_mg$dx %in% dx,]
  p_mg[[dx]] <- ggplot(dtss, aes(x = group, y = abs_freqdiff_dx2ct, fill = label)) +
    geom_bar(stat="identity", position = "dodge", color = NA) +
    ylim(0,28) +
    labs(x = "", y = "Abs. Diff. in mC Peak %", 
         title = paste(dx, ": Microglia-specific heterochromatin gain or loss 
       (empirical P <= 0.001 by Permutation 10k )", sep ="" )) +
    geom_text(size = 5, hjust=0.5, 
              #aes(label = fdr_label) ) +
              #aes(label = showP) ) +
              aes(label = P) ) +
    scale_fill_manual(values = colors.LOH) + add
  
}

fname <- paste("uniqueMC.ASC_MG.DXpeak_mCfreq.permutation.pdf", sep = ".")
pdf(file = file.path(outp, fname), width = 4, height = 2)
p_asc
p_mg
dev.off()  


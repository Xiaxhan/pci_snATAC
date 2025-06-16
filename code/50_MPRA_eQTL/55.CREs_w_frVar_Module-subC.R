##---- Identify modules based on Pattern of chromatin accessibility of (CREs contain SNPs) across subtypes/ subC-dx
##|a. Extract Peak Accessibility per pseudobulk (two versions: subC level and subC-dx level) from MACS callpeaks output
##|b. Depth normalization and log2-transofrm for each pseudobulk, then quantile normalize each peak across all pseudobulks
##|c. For subC (or subC-dx), calculate mean normalized peak score per peak across its pseudobulks
##|d. Subset peaks for CREs, Plot heatmap and show hiearchical clustering, cut dendrodam to define modules

library(data.table)
library(dplyr)
library(rhdf5)
library(parallel)
library(preprocessCore)
library(ComplexHeatmap)
library("RColorBrewer")
library("VennDiagram")
library(GenomicRanges)
library(factoextra)
outp <- "QTL_fsVars/out"

projmeta <- readRDS("projATAC_50subC_meta_584904.rds")
finalSubCs <- sort(unique(projmeta$subClusters))
asts <- colorRampPalette(c("darkred", "wheat"))(length(finalSubCs[finalSubCs %like% "ast"]))
mgs <- colorRampPalette(c("#F47D2B", "lightsalmon"))(length(finalSubCs[finalSubCs %like% "mg"]))
neus <- colorRampPalette(c("darkgreen","lightgreen"))(length(finalSubCs[finalSubCs %like% "neu"]))
odcs <- colorRampPalette(c("darkgoldenrod4", "lightgoldenrod1"))(length(finalSubCs[finalSubCs %like% "odc"]))
opcs <- colorRampPalette(c("cyan2", "lightcyan"))(length(finalSubCs[finalSubCs %like% "opc"]))
colors.subCs <- c(asts, mgs, neus, odcs, opcs)
names(colors.subCs) <- finalSubCs

colors.mg <-  ArchR::ArchRPalettes$stallion[1:10]
names(colors.mg) <- finalSubCs[finalSubCs %like% "mg"]

colors.dx <- ArchR::ArchRPalettes$circus[1:4]
names(colors.dx) <- c("Control","AD","bvFTD","PSP_S")

# ------------------Load and Make Objects
## 1.1 Load snp and peak table
mpra <- readRDS(file.path(outp, "obj.mpra_w_LDblock.rds"))
colnames(mpra)[4] <- "mpra_variant_id"
mpra <- mpra[!is.na(mpra$start_hg38), ]
#dt_peaks <- readRDS(file.path(outp, "obj.peaks_DynInv_max_mlogP_x_abslog2fc.rds"))

## 1.2 Load Peak annotaton "pciATAC_peakset"
source("pci_snATAC_peakSet_w_Enh.R")
pciATAC_peakset$peakloci <- paste(pciATAC_peakset@seqnames, 
                                  pciATAC_peakset@ranges@start,
                                  pciATAC_peakset@ranges@start + pciATAC_peakset@ranges@width - 1,
                                  sep = "-") 
dt_peakinfo <- values(pciATAC_peakset[,c("peakloci","peakType","nearestGene")])
rownames(dt_peakinfo) <- dt_peakinfo$peakloci

## 1.3 Load Peak Matrix  "QNmat_Pseudobulk_x_PE_subC" & "QNmat_Pseudobulk_x_PE_subCdx"
#   was generated in "/home/xhan/atac/archrproject0318/CRE/analyze/CREmodule_psedubulk.R" or "CREmodule_subCdx_pseudobulk.R"
projdir <- "out_rds_csv_projrmSubset_subPeak"

QNmat_Pseudobulk_x_PE_subC <- readRDS(file.path(projdir, "QuantileNorm_PseudobulkSubC_x_allpeaks.rds")) # 267 x 924225
QNmat_Pseudobulk_x_PE_subCdx <- readRDS(file.path(projdir, "QuantileNorm_PseudobulkSubCdx_x_allpeaks.rds")) # 823 x 924225
dim(QNmat_Pseudobulk_x_PE_subC)
dim(QNmat_Pseudobulk_x_PE_subCdx)

## 1.4 Extract snp-related peaks "dt_Peaks_w_snps"
allpeaks <- colnames(QNmat_Pseudobulk_x_PE_subC)

tmpV <- strsplit(allpeaks, split = "-") %>% unlist()
gr_allpeaks <- GRanges(
  seqnames=tmpV[seq(1, length(tmpV)-2, by = 3)],
  IRanges(start=as.numeric(tmpV[seq(2, length(tmpV)-1, by = 3)]),
          end=as.numeric(tmpV[seq(3, length(tmpV), by = 3)])
  )
)
values(gr_allpeaks) <- dt_peakinfo[allpeaks,]

used_snp <-  mpra[mpra$pval < 0.05,]
gr_snp <- GRanges(
  seqnames=paste("chr", used_snp$chr_hg38, sep = ""),
  IRanges(start=used_snp$start_hg38,
          end=used_snp$start_hg38+1)
)
values(gr_snp) <- used_snp

o <- findOverlaps(gr_snp, gr_allpeaks)
map_snp <- gr_snp[queryHits(o),]
map_peak <- gr_allpeaks[subjectHits(o),]

dt_Peaks_w_snps <- cbind(values(map_snp), values(map_peak))

v_CREs_w_snps <- unique(dt_Peaks_w_snps[dt_Peaks_w_snps$peakType %in% c("PE","Promoter"),]$peakloci) 

## 1.5 Load dx vs. ct max differential accessibility: 
# -log(Pval)*logFC        for marker peaks in per condition (keep up and down!)
# -log(Pval)*abs(logFC)   show the max changes across 4 conditions (ignore up/down)
diffscore_allpeaks <- readRDS(file.path(outp, "All_DiffScore_50subCs.rds"))
diffscore_allpeaks$peakloci <- gsub(":","-",diffscore_allpeaks$loci)

cre_diffscore <- diffscore_allpeaks[diffscore_allpeaks$peakloci %in% v_CREs_w_snps,]



## 1.6 load gwas snps
GWASsnp <-read.csv("/geschwindlabshares/RexachGroup/Xia_Data/resources/GWAS/FTD_META_bvFTDpval_updated_rsID_w_hg38_loci.txt", row.names = 1)
GWASsnp_wloci <- GWASsnp[!is.na(GWASsnp$hg38_loci), ]
dim(GWASsnp_wloci) #6,004,026

#---------1.7 Use Matrix for CREs x subC (or subCdx)
# pseudobulks mg-subC
if(T){
  mg_sps <- rownames(QNmat_Pseudobulk_x_PE_subC)[  rownames(QNmat_Pseudobulk_x_PE_subC) %like% "mg"]
  peakMat_subC <- QNmat_Pseudobulk_x_PE_subC[mg_sps,v_CREs_w_snps]
  dim(peakMat_subC) # 49 pseudobulks x 623 CREs
}

# pseudobulks mg-subCdx
if(T){
  mg_sps <- rownames(QNmat_Pseudobulk_x_PE_subCdx)[  rownames(QNmat_Pseudobulk_x_PE_subCdx) %like% "mg"]
  peakMat_subCdx <- QNmat_Pseudobulk_x_PE_subCdx[mg_sps,v_CREs_w_snps]
  dim(peakMat_subCdx) # 153 pseudobulks x 623 CREs
}

# mean per subC
avg_peakMat_subC <- NULL
for (subC in finalSubCs[finalSubCs %like% "mg"]){
  keepRows <- rownames(peakMat_subC)[  rownames(peakMat_subC) %like% subC ]
  pm <- peakMat_subC[keepRows, ]
  pm_mean <- colMeans(pm) 
  mat <- matrix(pm_mean, nrow = 1)
  colnames(mat) <- colnames(pm) 
  rownames(mat) <- subC  
  avg_peakMat_subC <- rbind(avg_peakMat_subC, mat)
}

# mean per subC-dx
avg_peakMat_subCdx <- NULL
for (subC in finalSubCs[finalSubCs %like% "mg"]){
  for (dx in c("Control","AD","bvFTD","PSP_S")){
    subC_dx <- paste(subC, dx, sep = ".")
    keepRows <- rownames(peakMat_subCdx)[  rownames(peakMat_subCdx) %like% subC_dx ]
    pm <- peakMat_subCdx[keepRows, ]
    pm_mean <- colMeans(pm) 
    mat <- matrix(pm_mean, nrow = 1)
    colnames(mat) <- colnames(pm) 
    rownames(mat) <- paste(subC, dx, sep = ".")
    avg_peakMat_subCdx <- rbind(avg_peakMat_subCdx, mat)
  }
}

# -----------------2. subC - Heatmap of Modules
useMat <- avg_peakMat_subC
matK_peakXsp_scaled <- t( scale(useMat) )

if(T){
  top_annotation_subc_only <- ComplexHeatmap::HeatmapAnnotation(
    subcluster = names(colors.mg),
    col = list(subcluster = colors.mg), #green
    which = "column")
  
  p_scale <- ComplexHeatmap::Heatmap(
    matK_peakXsp_scaled,
    #  col = colors,
    column_title = "623 sCREs_with_SNPs x 10 subCs, quantile norm, scale by peaks",
    na_col = "gray",
    column_title_rot = 0, # 0, 90, 270
    border = F, use_raster = T,
    show_row_names=F,show_column_names = T,
    cluster_columns = T, show_column_dend = T,
    cluster_rows = T,show_row_dend = T,
    column_names_side = "top",
    show_heatmap_legend = T,
    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
    #right_annotation = row_annotation,
    top_annotation = top_annotation_subc_only
  )
  p_scale <- draw(p_scale, heatmap_legend_side="bottom", annotation_legend_side="right",
                   legend_grouping = "original")
  
}

# Extract the dendrogram and cut the tree to define  modules
dendrogram_row <- p_scale@ht_list[[1]]@row_dend_list[[1]]
ordered_peaks <- p_scale@ht_list[[1]]@row_order
# convert the dendrogram into an hclust object
hclust_obj <- as.hclust(dendrogram_row)
# mat with peaks in heatmap-clustering order
ordered_mat <-  matK_peakXsp_scaled[ordered_peaks,]

# check 
pdf(file.path(outp,"dendrogram.Module.peak_w_frVar_x_subC.pdf"), width = 10, height = 8)
plot(dendrogram_row)
dev.off()

pdf(file.path(outp,"tt.pdf"), width = 10, height = 8)
ComplexHeatmap::Heatmap(ordered_mat, na_col = "gray",  column_names_side = "top",column_title_rot = 0, border = F, use_raster = T, show_row_names=F,show_column_names = T,cluster_columns = T, show_column_dend = T,cluster_rows = F,show_row_dend = F, show_heatmap_legend = T)
dev.off()

### Try for N Modules
cut_dendrom_plot <- function(tree_height){
  #row_clusters <- cutree(hclust_obj, k = 5)  # Cut by nb of clusters
  # row_clusters <- cutree(hclust_obj, h = 5.5) 
  row_clusters <- cutree(hclust_obj, h = tree_height) # Cut by tree height

  # Create row annotation (ensure this matches the new row order)
  modules <- unique(row_clusters)
  colors.modules <-  ArchR::ArchRPalettes$bear[1:length(modules)]
  names(colors.modules) <- modules
  
  row_annotation <- rowAnnotation(
    Module = anno_simple(row_clusters, col = colors.modules)
   # Text = anno_text(paste("Module", row_clusters), gp = gpar(fontsize = 10, col = "black"))  # Cluster labels as text
  )
  
  # Column annotation
  top_annotation_subc_only <- ComplexHeatmap::HeatmapAnnotation(
    subcluster = names(colors.mg),
    col = list(subcluster = colors.mg), #green
    which = "column")
  
  # color for values
  s_colors <- ArchR::ArchRPalettes$solarExtra
  breaks <- seq(min(matK_peakXsp_scaled), max(matK_peakXsp_scaled), length.out = length(s_colors))
  colors <- circlize::colorRamp2(breaks, s_colors)
  
  p_scale2 <- ComplexHeatmap::Heatmap(
    matK_peakXsp_scaled,
    #  col = colors,
    column_title = paste("623 sCREs_with_SNPs x 10 subCs, quantile norm, scale by peaks\n",
                         "#Cut row dendrogram height = ", tree_height, "\n",
                         "#Modules = ", length(modules), sep = ""),
    na_col = "gray",
    column_title_rot = 0, # 0, 90, 270
    border = F, use_raster = T,
    show_row_names=F,show_column_names = T,
    cluster_columns = T, show_column_dend = T,
    cluster_rows = T,show_row_dend = T,
    column_names_side = "top",
    show_heatmap_legend = T,
    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
    right_annotation = row_annotation,
    split = row_clusters, 
    top_annotation = top_annotation_subc_only

  )
  p_scale2 <- draw(p_scale2, heatmap_legend_side="bottom", annotation_legend_side="right",
                  legend_grouping = "original")
  
  #return(p_scale2)
}

pdf(file.path(outp,"heat.Module.peak_w_frVar_x_subC.pdf"), width = 10, height = 10)
p_scale
cut_dendrom_plot(tree_height=5)
cut_dendrom_plot(tree_height=5.2)
cut_dendrom_plot(tree_height=5.3)
cut_dendrom_plot(tree_height=5.4)
cut_dendrom_plot(tree_height=5.5)
cut_dendrom_plot(tree_height=5.6)
cut_dendrom_plot(tree_height=5.7)
dev.off()

# Cut tree and Plot heatmap with the cutting result
cut_dendrom_plot_splitRow <- function(tree_height){
  row_peaks <- rownames(ordered_mat)
  
  row_clusters <- cutree(hclust_obj, h = tree_height) # Cut by tree height
  row_clusters <- row_clusters[row_peaks]
  
  peak_to_module <- paste("M",row_clusters,sep = "")
  peak_to_module <- factor(peak_to_module, levels=unique(peak_to_module))
  
  # Create row annotation (ensure this matches the new row order)
  modules <- unique(row_clusters)
  colors.modules <-  ArchR::ArchRPalettes$bear[1:length(modules)]
  names(colors.modules) <- modules
  
  row_annotation <- rowAnnotation(
    Module = anno_simple(row_clusters, col = colors.modules)
  )
  
  # Column annotation
  top_annotation_subc_only <- ComplexHeatmap::HeatmapAnnotation(
    subcluster = names(colors.mg),
    col = list(subcluster = colors.mg), #green
    which = "column")
  
  # color for values
  s_colors <- ArchR::ArchRPalettes$solarExtra
  breaks <- seq(min(matK_peakXsp_scaled), max(matK_peakXsp_scaled), length.out = length(s_colors))
  colors <- circlize::colorRamp2(breaks, s_colors)
  
  p_scale2 <- ComplexHeatmap::Heatmap(
    ordered_mat,
    #  col = colors,
    column_title = paste("623 sCREs_with_SNPs x 10 subCs, quantile norm, scale by peaks\n",
                         "#Cut row dendrogram height = ", tree_height, "\n",
                         "#Modules = ", length(modules), sep = ""),
    na_col = "gray",
    column_title_rot = 0, # 0, 90, 270
    border = F, use_raster = T,
    show_row_names=F,show_column_names = T,
    cluster_columns = T, show_column_dend = T,
    cluster_rows = F,show_row_dend = F,
    column_names_side = "top",
    show_heatmap_legend = T,
    heatmap_legend_param = list(legend_direction = "horizontal", legend_width = unit(6, "cm")),
    right_annotation = row_annotation,
    row_split = peak_to_module, 
    top_annotation = top_annotation_subc_only
    
  )
  p_scale2 <- draw(p_scale2, heatmap_legend_side="bottom", annotation_legend_side="right",
                   legend_grouping = "original")
  
  #return(p_scale2)
}

pdf(file.path(outp,"heat.Module_SplitRow.peak_w_frVar_x_subC.pdf"), width = 10, height = 20)
p_scale
cut_dendrom_plot_splitRow(tree_height=5)
cut_dendrom_plot_splitRow(tree_height=5.2)
cut_dendrom_plot_splitRow(tree_height=5.3)
cut_dendrom_plot_splitRow(tree_height=5.4)
cut_dendrom_plot_splitRow(tree_height=5.5)
cut_dendrom_plot_splitRow(tree_height=5.6)
cut_dendrom_plot_splitRow(tree_height=5.7)
dev.off()


# -----------------3. Pick the cutting_height (5.2)
modules <- cutree(hclust_obj, h = 5.2) # Cut by tree height
df_module <- data.frame(modules)
colnames(df_module) <- "my_modules"

v_modules <- paste("module", modules, sep = "")
names(v_modules) <- names(modules)
table(v_modules) %>% data.frame() # save

cre_diffscore$module <- v_modules[cre_diffscore$peakloci]

# the order match heatmap with row tree (rowCluster = T)
p_pick <- cut_dendrom_plot(tree_height=5.2)
pick_ordered_peaks <- p_pick@ht_list[[1]]@row_order
pick_ordered_m <- p_pick@ht_list[[1]]@column_order 
ordered_peaks_module <-  rownames( matK_peakXsp_scaled[pick_ordered_peaks,] )
ordered_subc_module <-  colnames( matK_peakXsp_scaled[,pick_ordered_m] )
df_module[ordered_peaks_module,]

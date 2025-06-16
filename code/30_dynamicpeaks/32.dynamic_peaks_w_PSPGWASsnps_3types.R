#Generate a table of dynamic CREs associated with GWAS SNPs, 
#including SNP p-values and effect sizes (betas), using three criteria: 
#(1) CREs that directly contain GWAS SNPs, (2) CREs located within 2 kb of GWAS SNPs, 
# and (3) CREs that contain SNPs in linkage disequilibrium (LD) with GWAS SNPs

library(GenomicRanges)
library(data.table)
library(dplyr)
outp <- "out_rds_csv_projrmSubset_subPeak"
majorCs <- c("ast","mg","EX","IN","odc","opc")

projdir <- "outdir_projrmSubset_subPeak"
dynamicP_per_subc <- readRDS(file.path(projdir, "out_rds_csv_projrmSubset_subPeak","dynamicP_per_subc.rds") )

Dxs <- c("Control","AD","bvFTD","PSP_S")
keepCols <- c("peakloci",Dxs,"subC_dynP","TEclasses","nearestGene","PE_G","peakType","peakTypeRaw")
values(dynamicP_per_subc) <- values(dynamicP_per_subc)[,keepCols]

# PSP gwas
GWAS <- read.csv("PSP_Hoglinger_updated_w_hg38_loci.txt", sep = ",", row.names = 1)
head(GWAS)
rownames(GWAS) <- GWAS$SNP

gr_SNPs <- GRanges(seqnames = gsub(":.*","",GWAS$hg38_loci),
                   ranges = IRanges(as.numeric(gsub(".*\\:","",GWAS$hg38_loci)),
                                    as.numeric(gsub(".*\\:","",GWAS$hg38_loci))),
                   rsID = GWAS$SNP,
                   hg38_snploci = GWAS$hg38_loci,
                   pval = GWAS$P)

#1. Find SNPs locate within peaks
o <- findOverlaps(dynamicP_per_subc, gr_SNPs)

dt_snp_in_dynP <- cbind(
  values(dynamicP_per_subc[queryHits(o)]),
  values(gr_SNPs[subjectHits(o)])
)
rownames(dt_snp_in_dynP) <- NULL
dt_snp_in_dynP <- data.frame(dt_snp_in_dynP)

dim(dt_snp_in_dynP)  

dt_snp_in_dynP$label <- "within_peak"

#2. Find SNPs locate in peaks +- 2kb (not peak body)
dynamicP_extend <- dynamicP_per_subc
extend_bp <- 2000
dynamicP_extend <-  GRanges(seqnames = dynamicP_extend@seqnames,
                            ranges = IRanges(dynamicP_extend@ranges@start - extend_bp, 
                                             dynamicP_extend@ranges@start + dynamicP_extend@ranges@width - 1 + extend_bp
                            )) 
values(dynamicP_extend) <- values(dynamicP_per_subc)[,c(1:11)]
o2 <- findOverlaps(dynamicP_extend, gr_SNPs)

dt_snp_in_dynP_2kb <- cbind(
  values(dynamicP_extend[queryHits(o2)]),
  values(gr_SNPs[subjectHits(o2)])
)
dt_snp_in_dynP_2kb <- data.frame(dt_snp_in_dynP_2kb)

dim(dt_snp_in_dynP_2kb)  

dim(dt_snp_in_dynP) 

# Filter out snps(2kb) already identified within peaks
dt_snp_in_dynP_2kb$raw_snp_in_2kb <- dt_snp_in_dynP_2kb$rsID

dt_snp_in_dynP_2kb$rsID <- ifelse(dt_snp_in_dynP_2kb$rsID %in% dt_snp_in_dynP$rsID, NA, "new")
dt_snp_in_dynP_2kb$rsID <- ifelse(dt_snp_in_dynP_2kb$rsID %in% "new", dt_snp_in_dynP_2kb$raw_snp_in_2kb, NA)
head(dt_snp_in_dynP_2kb)
dt_snp_in_dynP$raw_snp_in_2kb <- NA

#
dt_snp_in_dynP_2kb$label <- "within_peak_2kb"

head(dt_snp_in_dynP)
head(dt_snp_in_dynP_2kb)

dynP_w_snps_2groups<- rbind(dt_snp_in_dynP,
                            dt_snp_in_dynP_2kb[,colnames(dt_snp_in_dynP_2kb)])

#3. Get SNPs in LD with GWAS SNPs 
dynamicP_ldsnp <- read.csv("ldscore-input.PSPgwas-LDsnps_in_dynP.csv")
dynamicP_ldsnp$X <- NULL

#4.Merge 3 kinds of snps
dynamicP_ldsnp$label <- "GWAS_LDsnp_in_peak"

setdiff(colnames(dynP_w_snps_2groups), colnames(dynamicP_ldsnp))
colnames(dynamicP_ldsnp)[colnames(dynamicP_ldsnp) %in% "SNP_B_PSPgwas"] <- "rsID"
dynamicP_ldsnp$hg38_snploci <- paste(dynamicP_ldsnp$seqnames, dynamicP_ldsnp$start, sep = ":")
dynamicP_ldsnp$raw_snp_in_2kb <- NA
dynamicP_ldsnp$pval <- GWAS[dynamicP_ldsnp$rsID,]$P

setdiff(colnames(dynamicP_ldsnp), colnames(dynP_w_snps_2groups))
dynP_w_snps_2groups$SNP_A_LDinput <- NA
dynP_w_snps_2groups$R2 <- NA

setdiff(colnames(dynP_w_snps_2groups), colnames(dynamicP_ldsnp))

dynamicP_related_snps <- rbind(dynP_w_snps_2groups,
                               dynamicP_ldsnp[,colnames(dynP_w_snps_2groups)])

colnames(dynamicP_related_snps)[colnames(dynamicP_related_snps) %in% "rsID"] <- "GWASsnp_or_ldSNPB"

astC1 <- dynamicP_related_snps[dynamicP_related_snps$subC_dynP %in% "ast.C1",]
write.csv(astC1, file.path(outp, "astC1_dynamicP-PSP-GWASsnps_3ways.txt"))

# save genes only for enrichment
# need R2 > 0.8
dt_genes <- dynamicP_related_snps[!(dynamicP_related_snps$label %in% "GWAS_LDsnp_in_peak" &
                                      dynamicP_related_snps$R2 <= 0.8), ]
dt_genes$GWASsnp_or_ldSNPB <- NULL
dt_genes$SNP_A_LDinput <- NULL
dt_genes$raw_snp_in_2kb <- NULL
dt_genes$R2 <- NULL

dt_genes$hg38_snploci <- NULL
dt_genes$pval_bvFTDgwas <- NULL
dt_genes$pval_metaFTDgwas <- NULL

# can filter by Dx & subC_dyn & label for gene sets
dt_genes <- unique(dt_genes)

write.csv(dt_genes, file.path(outp, "dynamicP_related_PSP-GWASsnps_3ways_Genes_R2_0.8.txt"))


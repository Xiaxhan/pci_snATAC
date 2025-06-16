library(memes)
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
PATHDB<- "/geschwindlabshares/RexachGroup/SharedData/motif_databases"
motifdb <- file.path(PATHDB, "CIS-BP_2.00/Homo_sapiens.meme")
memeBin <- "/home/xhan/apps/meme-5.4.1/bin"

run_TFenrich <- function(gr, label){
  seqs <- get_sequence(gr, genome)

  dir_n <- file.path(outp,  paste("memes",label,sep = "_"))
  if (!dir.exists(dir_n)){dir.create(dir_n)}
  runAme(input = seqs,
         control = "shuffle", # seq or "shuffle"
         outdir = dir_n,
         method = "fisher",
         database = motifdb,
         meme_path = memeBin)
}

run_TFenrich(gr_m4, "top_10pcf_mgC4_x_ftdGWAS")
run_TFenrich(gr_a1, "top_10pcf_astC1_x_pspGWAS")


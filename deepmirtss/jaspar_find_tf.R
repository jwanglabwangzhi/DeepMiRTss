args<-commandArgs(TRUE)
fasta_file_name <- args

library(TFBSTools)
library(JASPAR2018)
library(Biostrings)

opts <- list()
opts[["species"]] <- "9606"
#opts[["type"]] <- "ChIP-seq"  不选就可以选中所有type
pfms <- getMatrixSet(JASPAR2018, opts)
pwms <- toPWM(pfms)
fasta_read<- readDNAStringSet(fasta_file_name)
sitesetList <- searchSeq(pwms, fasta_read, strand="+", min.score="80%")
output_file_name <-"tf_jaspar.txt"
sink(output_file_name)
write.csv(writeGFF3(sitesetList), file=output_file_name)

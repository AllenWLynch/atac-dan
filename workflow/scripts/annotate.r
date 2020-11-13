#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) == 3){
  bedfile <- args[1]
  species <- args[2]
  output <- args[3]
} else {
  stop("Requires [bedfile] [species] [output]", call.=FALSE)
}

library(ChIPseeker)

if (reference == 'hg38'){
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
} else if (reference == 'mm10'){
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  txdb = TxDb.Mmusculus.UCSC.mm10.knownGene
} else {
  stop("Species must be either hg38 or mm10", call.=False)
}

peaks <- readPeakFile(bedfile)

annotations <- annotatePeak(peaks, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)

write.table(as.data.frame(annotations@anno), file = output, quote=FALSE, sep = "\t", row.names=FALSE)
#set data source

library(seqinr)
library(dplyr)
library(signalHsmm)

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

signalHsmm1987 <- train_hsmm(read_uniprot(paste0(pathway, "sp1950_1989.txt"), ft_names = "signal"),
                              aaaggregation)

plas_pos <- read_uniprot("plasmodium_pos.txt", euk = TRUE, what = "signal")
plas_neg <- read.fasta("plas_neg.fasta", seqtype = "AA")

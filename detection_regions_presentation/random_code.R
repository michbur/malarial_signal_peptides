library(ggplot2)
library(pbapply)
library(seqinr)
library(dplyr)
library(reshape2)
library(signalHsmm)
library(biogram)

seq50_10 <- read_uniprot("./training_data/sp1950_2010.txt", ft_names = "signal")

longer_seqs <- seq50_10[lengths(seq50_10) > 100]



mat_seq <- do.call(rbind, lapply(longer_seqs, function(single_seq)
  c(single_seq[(attr(single_seq, "signal")[2] + 1):100], rep(NA, attr(single_seq, "signal")[2]))))

sig_seq <- do.call(rbind, lapply(longer_seqs, function(single_seq)
  c(single_seq[2L:attr(single_seq, "signal")[2]], rep(NA, 101 - attr(single_seq, "signal")[2]))))

seqs <- rbind(sig_seq, mat_seq)
tar <- sapply(c(1, 0), function(i) rep(i, nrow(sig_seq)))
ngrams <- count_multigrams(ns = c(1, rep(2, 4), rep(3, 3)),
                           ds = list(0, 0, 1, 2, 3, c(0, 0), c(0, 1), c(1, 0)), seqs, u = a()[-1])

ngrams <- as.matrix(ngrams) > 0
storage.mode(ngrams) <- "integer"
tf_res <- test_features(tar, ngrams)
imps <- cut(tf_res, breaks = c(0, 0.00001, 1))[[1]]

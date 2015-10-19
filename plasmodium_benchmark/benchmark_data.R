# library(seqinr)
# dat <- read.fasta("benchmark_plas_data.fasta")
# #truncated data for signalP 3.0
# short_dat <- lapply(dat, function(i) i[1L:200])
# write.fasta(short_dat, names = names(dat), file.out = "benchmark_plas_data_tr.fasta")

library(hmeasure)
real_labels <- c(rep(1, 102), rep(0, 358))

HMeasure(real_labels, read.csv2("./plasmodium_benchmark/benchmark_plas_other.csv"))[["metrics"]]

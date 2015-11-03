# library(seqinr)
# dat <- read.fasta("benchmark_plas_data.fasta")
# #truncated data for signalP 3.0
# short_dat <- lapply(dat, function(i) i[1L:ifelse(length(i) < 200, length(i), 200)])
# write.fasta(short_dat, names = names(dat), file.out = "benchmark_plas_data_tr.fasta")

# library(hmeasure)
# real_labels <- c(rep(1, 102), rep(0, 358))
# 
# HMeasure(real_labels, read.csv2("./plasmodium_benchmark/benchmark_plas_other.csv"))[["metrics"]]

library(seqinr)
library(hmeasure)
library(signalHsmm)
source("./functions/benchmark_functions.R")



#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[19500000 TO 19870000] AND reviewed:yes
signalHsmm1987 <- train_hsmm(read_uniprot("./signalHsmm_1987_data/sp1950_1987.txt", ft_names = "signal"),
                             aaaggregation)


# BENCHMARK - PLASMODIUM -------------------------------------------------

#taxonomy:"Plasmodium [5820]" annotation:(type:signal evidence:manual) AND reviewed:yes
#111 proteins total, 102 proteins after purification
#taxonomy:"Plasmodium [5820]" NOT annotation:(type:signal evidence:manual) AND reviewed:yes
#361 proteins, 358 proteins after purification
real_labels <- c(rep(1, 102), rep(0, 358))



preds <- read_other_software("./plasmodium_benchmark_results")
signalHsmm1987_preds <- pred2df(predict(signalHsmm1987, read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta",
                                                                    seqtype = "AA")))[, "sp.probability"]
signalHsmm_preds <- pred2df(run_signalHsmm(read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta",
                                                      seqtype = "AA")))[, "sp.probability"]

calc_metrics(real_labels, data.frame(preds, signalHsmm = signalHsmm_preds, signalHsmm1987 = signalHsmm1987_preds), 0.005)

# library(ggplot2)
# dat <- data.frame(real = real_labels, pred = signalHsmm_preds)
# 
# ggplot(dat, aes(x = as.factor(real), y = pred)) +
#   geom_boxplot()
# 
# 
# library(dplyr)
# group_by(dat, real) %>% summarise(median(pred))

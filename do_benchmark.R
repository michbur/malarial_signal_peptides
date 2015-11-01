source("./functions/benchmark_functions.R")


# library(seqinr)
# dat <- read.fasta("benchmark_plas_data.fasta")
# #truncated data for signalP 3.0
# short_dat <- lapply(dat, function(i) i[1L:ifelse(length(i) < 200, length(i), 200)])
# write.fasta(short_dat, names = names(dat), file.out = "benchmark_plas_data_tr.fasta")

# library(hmeasure)
# real_labels <- c(rep(1, 102), rep(0, 358))
# 
# HMeasure(real_labels, read.csv2("./plasmodium_benchmark/benchmark_plas_other.csv"))[["metrics"]]




#taxonomy:"Plasmodium [5820]" annotation:(type:signal evidence:manual) AND reviewed:yes
#111 proteins
#taxonomy:"Plasmodium [5820]" NOT annotation:(type:signal evidence:manual) AND reviewed:yes
#361 proteins

read_other_software("./plasmodium_benchmark_results")

real_labels <- c(rep(1, 102), rep(0, 358))

HMeasure(real_labels, preds)

signal.hsmm1987_preds <- pred2df(predict(signal.hsmm1987, c(plas_pos, plas_neg)))
signal.hsmm_preds <- pred2df(run_signalHsmm(c(plas_pos, plas_neg)))



TP <- as.numeric(metrics[["TP"]])
FP <- as.numeric(metrics[["FP"]])
TN <- as.numeric(metrics[["TN"]])
FN <- as.numeric(metrics[["FN"]])

MCC = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))


other_software <- t(do.call(rbind, lapply(list(signalPnotm = read_signalp41("./benchmark/signaP41notm.txt"), 
                                      signalPtm = read_signalp41("./benchmark/signaP41tm.txt"), 
                                      predsi = read_predsi("./benchmark/predsi.txt"),
                                      phobius = read_phobius("./benchmark/phobius.txt")), function(predictor)
                                        predictor[["sp.probability"]])))
#write.csv2(other_software, file = "other_soft.csv")

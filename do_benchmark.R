library(seqinr)
library(hmeasure)
library(signalHsmm)
source("./functions/benchmark_functions.R")



# signalHsmm1986 and signalHsmm2010 -------------------------------------

#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[19500000 TO 19870000] AND reviewed:yes
#354 proteins, 335 after purification
signalHsmm1986 <- train_hsmm(read_uniprot("./training_data/sp1950_1987.txt", ft_names = "signal"),
                             aaaggregation)
#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[19500000 TO 19870000] AND reviewed:yes
#2372 proteins, 2313 after purification
signalHsmm2010 <- train_hsmm(read_uniprot("./training_data/sp1950_2010.txt", ft_names = "signal"),
                             aaaggregation)

# BENCHMARK - PLASMODIUM -------------------------------------------------

#taxonomy:"Plasmodium [5820]" annotation:(type:signal evidence:manual) AND reviewed:yes
#111 proteins total, 102 proteins after purification
#taxonomy:"Plasmodium [5820]" NOT annotation:(type:signal evidence:manual) AND reviewed:yes
#361 proteins, 358 proteins after purification

calc_metrics(c(rep(1, 102), rep(0, 358)), 
             data.frame(read_other_software("./plasmodium_benchmark_results"), 
                        signalHsmm2010 = pred2df(predict(signalHsmm2010, 
                                                         read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta",
                                                                    seqtype = "AA")))[, "sp.probability"], 
                        signalHsmm1987 = pred2df(predict(signalHsmm1986, 
                                                         read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta",
                                                                    seqtype = "AA")))[, "sp.probability"]), 0.5)

# BENCHMARK - all -------------------------------------------------

calc_metrics(c(rep(1, 214), rep(0, 214)), 
             data.frame(read_other_software("./benchmark_results"), 
                        signalHsmm2010 = pred2df(predict(signalHsmm2010, 
                                                         read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta",
                                                                    seqtype = "AA")))[, "sp.probability"], 
                        signalHsmm1987 = pred2df(predict(signalHsmm1986, 
                                                         read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta",
                                                                    seqtype = "AA")))[, "sp.probability"]), 0.5)

# library(ggplot2)
# dat <- data.frame(real = real_labels, pred = signalHsmm_preds)
# 
# ggplot(dat, aes(x = as.factor(real), y = pred)) +
#   geom_boxplot()
# 
# 
# library(dplyr)
# group_by(dat, real) %>% summarise(median(pred))

# source of proteins in signalHsmm1986
# prot86 <- read_uniprot("./training_data/sp1950_1987.txt", ft_names = "signal")
# table(unname(sapply(prot86, function(i) attr(i, "OS"))))
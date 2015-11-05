#set data source

library(seqinr)
library(dplyr)
library(signalHsmm)
library(hmeasure)
library(pbapply)
library(reshape2)

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

source("./functions/plot_tools.R")
source("./functions/benchmark_functions.R")

# CROSS-VALIDATION

load(paste0(pathway, "fold_res_df.RData"))




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

metrics_plas <- calc_metrics(c(rep(1, 102), rep(0, 358)), 
                             data.frame(read_other_software("./plasmodium_benchmark_results"), 
                                        signalHsmm2010 = pred2df(predict(signalHsmm2010, 
                                                                         read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta",
                                                                                    seqtype = "AA")))[, "sp.probability"], 
                                        signalHsmm1987 = pred2df(predict(signalHsmm1986, 
                                                                         read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta",
                                                                                    seqtype = "AA")))[, "sp.probability"]), 0.005)

# BENCHMARK - ALL -------------------------------------------------

metrics_all <- calc_metrics(c(rep(1, 214), rep(0, 214)), 
                            data.frame(read_other_software("./benchmark_results"), 
                                       signalHsmm2010 = pred2df(predict(signalHsmm2010, 
                                                                        read.fasta("./benchmark_data/benchmark_data.fasta",
                                                                                   seqtype = "AA")))[, "sp.probability"], 
                                       signalHsmm1987 = pred2df(predict(signalHsmm1986, 
                                                                        read.fasta("./benchmark_data/benchmark_data.fasta",
                                                                                   seqtype = "AA")))[, "sp.probability"]), 0.005)
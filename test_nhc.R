require(XML)
require(seqinr)
require(fitdistrplus)
require(dplyr)
require(signalHsmm)
require(hmeasure)
require(pbapply)
require(reshape2)
require(hmeasure)
require(xtable)
require(biogram)
require(ggplot2)
require(grid)
require(gridExtra)

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

source("./functions/plot_tools.R")
source("./functions/reglen_plot.R")
source("./functions/cv_analysis.R")
source("./functions/enc_region.R")
source("./functions/benchmark_functions.R")
source("./functions/cdhit.R")
source("./functions/train_signalHsmms.R")
source("./functions/signalHsmm_kmer.R")

# absurdly simple
find_nhc2 <- function(protein, signal = NULL) {
  protein <- toupper(protein)
  if (is.null(signal)) 
    signal <- attr(protein, "signal")
  
  sig <- protein[signal[1]:signal[2]]
  
  start_c <- signal[2] - 4
  
  start_n <- 1
  start_h <- 2
  noh <- 0
  
  while(noh < 4 && start_h < start_c) {
    start_h <- start_h + 1
    noh <- noh + ifelse(sig[start_h] %in% c("F", "I", "L", "V", "W"), 1, 0)
    noh <- ifelse(noh < 1, 0, noh)
  }
  
  start_h <- start_h - 4
  
  c(start_n = 1, start_h = start_h, start_c = start_c, cs = signal[2])
}

find_nhc3 <- function(protein, signal = NULL) {
  protein <- toupper(protein)
  if (is.null(signal)) 
    signal <- attr(protein, "signal")
  
  sig <- protein[signal[1]:signal[2]]
  
  start_c <- signal[2] - 4
  
  start_h <- start_c - 6
  #if statement to prevent negative values
  if (start_h > 1) {
    #nonh number of nonhydrophobic residues
    nonh <- 0
    #noc number of charged
    noc <- 0
    while(nonh < 3 && noc == 0 && start_h > 1) {
      start_h <- start_h - 1
      nonh <- nonh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), -1, 1)
      nonh <- ifelse(nonh < 0, 0, nonh)
      noc <- ifelse(sig[start_h] %in% c("R", "H", "K", "D", "E"), 1, 0)
    }
  } else {
    start_h <- 1
  }
  
  prestart_c <- start_c - 1
  noh <- 0
  while(noh == 0 && start_h < prestart_c) {
    start_h <- start_h + 1
    noh <- noh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, 0)
  }
  #c(start_n = signal[1], start_h = start_h, start_c = start_c, cs = signal[2])
  c(start_n = 1, start_h = start_h, start_c = start_c, cs = signal[2])
}

# sapply(seq50_10, find_nhc2)
# find_nhc2(seq50_10[[2031]])

train_signalHsmms_nhc <- function(seqs) {
  #iterations without degeneration
  #aas <- tolower(a()[-1])
  aas <- a()[-1]
  names(aas) <- 1L:20
  
  # no degeneration
  signalHsmmNODEG <- train_hsmm(seqs, aas)
  signalHsmmNODEG2 <- train_hsmm(seqs, aas, region_fun = find_nhc2)
  signalHsmmNODEG3 <- train_hsmm(seqs, aas, region_fun = find_nhc3)
  
  list(signalHsmmNODEG = signalHsmmNODEG,
       signalHsmmNODEG2 = signalHsmmNODEG2,
       signalHsmmNODEG3 = signalHsmmNODEG3)
}


# removed single 7 aa long presequence
signalHsmms10 <- train_signalHsmms_nhc(seq50_10[-2031])
names(signalHsmms10) <- paste0(names(signalHsmms10), "_10")

signalHsmms87 <- train_signalHsmms_nhc(seq50_87)
names(signalHsmms87) <- paste0(names(signalHsmms87), "_87")

other_pred_plas <- read_other_software("./plasmodium_benchmark_results_NOHOM")
signalHsmm_pred_plas <- get_signalHsmm_preds(c(signalHsmms10, signalHsmms87),
                                             "./plasmodium_benchmark_data/benchmark_plas_data_NOHOM.fasta")

metrics_plas_NOHOM <- calc_metrics(c(rep(1, 51), rep(0, 211)), 
                                   data.frame(other_pred_plas[["prob"]], signalHsmm_pred_plas[["prob"]]), 
                                   0.5)
metrics_plas_NOHOM[c("AUC", "MCC")]

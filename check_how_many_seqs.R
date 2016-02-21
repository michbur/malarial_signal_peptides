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

check_how_many <- function(seqs) {
  #homology 50
  seq50f <- cdhit(seqs, thresh = 0.5, word_length = 2, only_signal = TRUE)
  
  #homology 90
  seq90f <- cdhit(seqs, thresh = 0.9, word_length = 5, only_signal = TRUE)
  
  c(seq50f = length(seq50f), seq90f = length(seq90f))
}

# data
#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[19500000 TO 19870000] AND reviewed:yes
#354 proteins, 335 after purification
seq50_87 <- read_uniprot("./training_data/sp1950_1987.txt", ft_names = "signal")

#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[19500000 TO 20100000] AND reviewed:yes
#2372 proteins, 2313 after purification
seq50_10 <- read_uniprot("./training_data/sp1950_2010.txt", ft_names = "signal")

check_how_many(seq50_10)

check_how_many(seq50_87)


#set data source

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
library(ggplot2)
library(grid)
library(gridExtra)

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"

source("./functions/plot_tools.R")
source("./functions/reglen_plot.R")
source("./functions/cv_analysis.R")
source("./functions/enc_region.R")
source("./functions/benchmark_functions.R")




# DENSITY OF LENGTH DISTRIBUTION (BOTH SIGNAL PEPTIDES AND THEIR REGIONS) --------------------------

reglen <- plot_reglen()

cairo_ps("./publication/figures/reglen.eps", width = 9, height = 5, onefile = FALSE)
grid.arrange(textGrob("A", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), reglen[["sp_len"]], 
             textGrob("B", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), reglen[["regions"]], 
             nrow = 2, ncol = 2, widths = c(0.05, 0.95))
dev.off()


# CROSS-VALIDATION -------------------------

# time-consuming part of the code silenced
# load("/Qnap/Publikacje/gcbsubmission_data/signalHsmm_cv.RData")
# rep_res <- perf_rep(fold_res, 0.05)
# save(rep_res, file = "./analysis_data/cv_results.RData")
load("./analysis_data/cv_results.RData")
cvplot <- plot_cvplot(create_cvplotdat(rep_res))

cairo_ps("./publication/figures/cvres.eps", width = 9, height = 5, onefile = FALSE)
print(cvplot[["plot"]])
dev.off()
cvplot[["cpt"]]
cat(cvplot[["xtab"]])

# THE BEST ENCODINGS -------------------------


enc_region <- create_enc_region(p1_dat = create_cvplotdat(rep_res))

cairo_ps("./publication/figures/enccomp.eps", width = 9, height = 8, onefile = FALSE)
grid.arrange(textGrob("A", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), enc_region[["prop_plot"]], 
             textGrob("B", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), enc_region[["freq_plot"]],
             nrow = 2, ncol = 2, widths = c(0.05, 0.95))
dev.off()

cat(enc_region[["best_sens"]])
cat(enc_region[["best_spec"]])

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
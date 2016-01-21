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

# data
#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[19500000 TO 19870000] AND reviewed:yes
#354 proteins, 335 after purification
seq50_87 <- read_uniprot("./training_data/sp1950_1987.txt", ft_names = "signal")

#taxonomy:"Eukaryota [2759]" annotation:(type:signal evidence:experimental) created:[19500000 TO 20100000] AND reviewed:yes
#2372 proteins, 2313 after purification
seq50_10 <- read_uniprot("./training_data/sp1950_2010.txt", ft_names = "signal")

signalHsmms10 <- train_signalHsmms(seq50_10)
names(signalHsmms10) <- paste0(names(signalHsmms10), "_10")

signalHsmms87 <- train_signalHsmms(seq50_87)
names(signalHsmms87) <- paste0(names(signalHsmms87), "_87")

# save(signalHsmms87, signalHsmms10, file = "./cache/signalHsmms.RData")
# load("./cache/signalHsmms.RData")

# BENCHMARK - PLASMODIUM HOMOLOGY REDUCED -------------------------------------------------

#taxonomy:"Plasmodium [5820]" annotation:(type:signal evidence:manual) AND reviewed:yes
#111 proteins total, 102 proteins after purification
#taxonomy:"Plasmodium [5820]" NOT annotation:(type:signal evidence:manual) AND reviewed:yes
#361 proteins, 358 proteins after purification

# all_seqs_plas <- read.fasta("./plasmodium_benchmark_data/benchmark_plas_data.fasta", seqtype = "AA")
# all_seqs_plasf <- cdhit(all_seqs_plas, thresh = 0.5, word_length = 2, only_signal = FALSE)
# write.fasta(lapply(all_seqs_plas[all_seqs_plasf], function(i) i[1L:150]), 
#             names = all_seqs_plasf, 
#             file.out = "./plasmodium_benchmark_data/benchmark_plas_data_NOHOM.fasta")

# et <- c(rep(1, 102), rep(0, 358))
# names(et) <- names(all_seqs_plas)
# etiquettes are commented out and written by hand in calc_metrics
# et[all_seqs_plasf]

metrics_plas_NOHOM <- calc_metrics(c(rep(1, 51), rep(0, 211)), 
                                   data.frame(read_other_software("./plasmodium_benchmark_results_NOHOM")[["prob"]], 
                                              #signalHsmmKmer_10 = signalHsmm_kmer(signalHsmms10[["signalHsmm_10"]], 
                                              #                                    read.fasta("./plasmodium_benchmark_data/benchmark_plas_data_NOHOM.fasta", seqtype = "AA")),
                                              get_signalHsmm_preds(c(signalHsmms10, signalHsmms87),
                                                                   "./plasmodium_benchmark_data/benchmark_plas_data_NOHOM.fasta")[["prob"]]
                                              ), 
                                   0.5)

write.csv(round(metrics_plas_NOHOM, 6), file = "./publication/supplements/S1_plasmodium_benchmark.csv")

pub_tab <- format_bench_table(metrics_plas_NOHOM, 
                              caption = "Comparison of Area Under the Curve, H-measure and Matthews Correlation Coefficient 
for different classifiers considering proteins belonging to \\textit{Plasmodiidae}.",
                              "tab:bench2010plas")
cat(pub_tab)



# BENCHMARK - ALL TAXONS NO HOMOLOGOUS -------------------------------------

# removal of homologs is commented out (time consuming)
# sp_seqs <- read.fasta("./benchmark_data/sp2010_2015.fasta", seqtype = "AA")
# nsp_seqs <- read.fasta("./benchmark_data/nsp2010_2015.fasta", seqtype = "AA")
# 
# sp_seqsf <- cdhit(sp_seqs, thresh = 0.5, word_length = 2, only_signal = FALSE)
# nsp_seqsf <- cdhit(nsp_seqs, thresh = 0.5, word_length = 2, only_signal = FALSE)

# set.seed(1)
# chosen_nsp <- nsp_seqsf[sample(1L:length(nsp_seqsf), length(sp_seqsf), replace = FALSE)]
# seq_NOHOM <- c(sp_seqs[sp_seqsf], nsp_seqs[chosen_nsp])
# 
# write.fasta(lapply(seq_NOHOM, function(i) i[1L:ifelse(length(i) > 150, 150, length(i))]), 
#             names = names(seq_NOHOM), 
#             file.out = "./benchmark_data/benchmark_data_NOHOM.fasta")

# metrics_all_NOHOM <- calc_metrics(c(rep(1, length(sp_seqsf)), rep(0, length(sp_seqsf))), 
metrics_all_NOHOM <- calc_metrics(c(rep(1, 127), rep(0, 127)), 
                                  data.frame(read_other_software("./benchmark_results_NOHOM")[["prob"]],
                                             #signalHsmmKmer_10 = signalHsmm_kmer(signalHsmms10[["signalHsmm_10"]], 
                                             #                                    read.fasta("./benchmark_data/benchmark_data_NOHOM.fasta", seqtype = "AA")),
                                             get_signalHsmm_preds(c(signalHsmms10, signalHsmms87),
                                                                  "./benchmark_data/benchmark_data_NOHOM.fasta")[["prob"]]
                                             ), 
                                  0.5)

write.csv(round(metrics_all_NOHOM, 6), file = "./publication/supplements/S2_general_benchmark.csv")

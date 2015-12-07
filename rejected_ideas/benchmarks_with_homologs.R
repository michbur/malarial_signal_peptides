# BENCHMARK - PLASMODIUM -------------------------------------------------

metrics_plas <- calc_metrics(c(rep(1, 102), rep(0, 358)), 
                             data.frame(read_other_software("./plasmodium_benchmark_results"), 
                                        get_signalHsmm_preds(c(signalHsmms10, signalHsmms87),
                                                             "./plasmodium_benchmark_data/benchmark_plas_data.fasta")), 0.005)

metrics_all <- calc_metrics(c(rep(1, 214), rep(0, 214)), 
                            data.frame(read_other_software("./benchmark_results"), 
                                       get_signalHsmm_preds(c(signalHsmms10, signalHsmms87),
                                                            "./benchmark_data/benchmark_data.fasta")), 0.005)
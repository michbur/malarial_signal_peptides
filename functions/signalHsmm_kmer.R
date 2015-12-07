signalHsmm_kmer <- function(signalHsmm_model, benchmark_data) {
  
  add_k_mer_state_easy <- function(signalHsmm_model, kMer1, pState1, nState1, pTrans1, d)
    add_k_mer_state(kMer1, signalHsmm_model[["pipar"]], signalHsmm_model[["tpmpar"]], 
                    signalHsmm_model[["od"]], signalHsmm_model[["params"]], 
                    pState1, nState1, pTrans1, d)
  
  maxSignal <- 45
  
  p4 <- 0.1
  p3 <- 0.25/(1-p4)
  p2 <- 0.35/(1-p4)/(1-p3)
  p1 <- 0.1/(1-p4)/(1-p3)/(1-p2)
  
  kMer1 <- list(3, 1, 2, 1, c(1,3))
  pState1 <- 3
  nState1 <- 4
  pTrans1 <- p1
  
  # -3|222
  kMer2 <- list(c(2,3),c(1,2), c(2,4))
  pState2 <- 3
  nState2 <- 4
  pTrans2 <- p2
  
  #-3|4(1,2,3,4)(2,4)
  kMer3 <- list(4, c(1,2,3), c(2,4))
  pState3 <- 3
  nState3 <- 4
  pTrans3 <- p3
  
  # kMer4
  kMer4 <- list(2, 1, c(1,4))
  pState4 <- 3
  nState4 <- 4
  pTrans4 <- p4
  
  lastState1 = 4 + 1 + length(kMer1)
  lastState2 = lastState1 + 1 + length(kMer2)
  lastState3 = lastState2 + 1 + length(kMer3)
  lastState4 = lastState3 + 1 + length(kMer4)
  
  parametersSet5 <- signalHsmms10[["signalHsmm_10"]] %>%
    add_k_mer_state_easy(kMer1, pState1, nState1, pTrans1, d = 3) %>%
    add_k_mer_state_easy(kMer2, pState2, nState2, pTrans2, d = 3) %>%
    add_k_mer_state_easy(kMer3, pState3, nState3, pTrans3, d = 3) %>%
    add_k_mer_state_easy(kMer4, pState4, nState4, pTrans4, d = 3)
  
  overall_probs_log <- signalHsmm_model$overall_probs_log
  maxSignal <- 45
  
  #registerDoMC(cores=4)
  n <- length(benchmark_data)
  results <- lapply(1L:n, function(i) {
    prot <- benchmark_data[[i]]
    deg_sample <- na.omit(as.numeric(degenerate(toupper(prot)[1L:min(maxSignal, length(prot))], aaaggregation)))
    
    viterbi_res <- duration_viterbi(deg_sample-1, signalHsmm_model$pipar, signalHsmm_model$tpmpar, 
                                    signalHsmm_model$od, signalHsmm_model$params)
    viterbi_path <- viterbi_res[["path"]]+1
    c_site <- ifelse(any(viterbi_path == 4), max(which(viterbi_path == 3)), length(deg_sample))
    prob.signal <- viterbi_res[["viterbi"]][c_site, viterbi_path[c_site]]
    prob.non <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site], 0)
    prob.total <- exp(prob.signal - prob.non)
    
    viterbi_res5 <- duration_viterbi(deg_sample-1, parametersSet5$pipar, parametersSet5$tpmpar, 
                                     parametersSet5$od, parametersSet5$params)
    viterbi_path5 <- viterbi_res5[["path"]]+1
    c_site5 <- ifelse(any(viterbi_path5 == 4), min(which(viterbi_path5 == 4))-1, length(deg_sample))
    c_site5 <- ifelse(any(viterbi_path5 == lastState1), min(which(viterbi_path5 == lastState1)), c_site5)
    c_site5 <- ifelse(any(viterbi_path5 == lastState2), min(which(viterbi_path5 == lastState2)), c_site5)
    c_site5 <- ifelse(any(viterbi_path5 == lastState3), min(which(viterbi_path5 == lastState3)), c_site5)
    c_site5 <- ifelse(any(viterbi_path5 == lastState4), min(which(viterbi_path5 == lastState4)), c_site5)
    prob.signal5 <- viterbi_res5[["viterbi"]][c_site5, viterbi_path5[c_site5]]
    prob.non5 <- Reduce(function(x, y) x + overall_probs_log[y], deg_sample[1L:c_site5], 0)
    prob.total5 <- exp(prob.signal5 - prob.non5)
    
    c(1 - 1/(1 + prob.total), 1 - 1/(1 + prob.total5), c_site, c_site5)
  })
  probs <- t(sapply(results, function (x) head(x, 2)))
  cut <- t(sapply(results, function (x) tail(x, 2)))
  probs[, 2]
}
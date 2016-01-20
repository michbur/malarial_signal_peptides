#' Train iterations of signalHsmm
#'
#' Trains all iterations of signalHsmm.
#'
#' @param seqs list of sequences.
#'
#' @return A list of signalHsmm predictors.

train_signalHsmms <- function(seqs) {
  #iterations without degeneration
  #aas <- tolower(a()[-1])
  aas <- a()[-1]
  names(aas) <- 1L:20
  
  # no degeneration
  signalHsmmNODEG <- train_hsmm(seqs, aas)
  signalHsmmNOCNODEG <- train_hsmm(seqs, aas, region_fun = find_nhc2)
  
  # degeneration
  signalHsmm <- train_hsmm(seqs, aaaggregation)
  signalHsmmNOC <- train_hsmm(seqs, aaaggregation, region_fun = find_nhc2)
  
  #homology 50
  seq50f <- cdhit(seqs, thresh = 0.5, word_length = 2, only_signal = TRUE)
  signalHsmmNOHOM50 <- train_hsmm(seqs[seq50f], aaaggregation)
  signalHsmmNOCNOHOM50 <- train_hsmm(seqs, aaaggregation, region_fun = find_nhc2)
  
  #homology 90
  seq90f <- cdhit(seqs, thresh = 0.9, word_length = 5, only_signal = TRUE)
  signalHsmmNOHOM90 <- train_hsmm(seqs[seq90f], aaaggregation)
  signalHsmmNOCNOHOM90 <- train_hsmm(seqs[seq90f], aaaggregation, region_fun = find_nhc2)
  
  list(signalHsmm = signalHsmm,
       signalHsmmNODEG = signalHsmmNODEG,
       signalHsmmNOHOM50 = signalHsmmNOHOM50,
       signalHsmmNOHOM90 = signalHsmmNOHOM90,
       signalHsmmNOC = signalHsmmNOC,
       signalHsmmNOCNODEG = signalHsmmNOCNODEG,
       signalHsmmNOCNOHOM50 = signalHsmmNOCNOHOM50,
       signalHsmmNOCNOHOM90 = signalHsmmNOCNOHOM90)
}
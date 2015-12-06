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

  # no degeneration
  signalHsmm <- train_hsmm(seqs, aaaggregation)
  
  #homology 50
  seq50f <- cdhit(seqs, thresh = 0.5, word_length = 2, only_signal = TRUE)
  signalHsmmNOHOM50 <- train_hsmm(seqs[seq50f], aaaggregation)
  
  #homology 90
  seq90f <- cdhit(seq50_10, thresh = 0.9, word_length = 5, only_signal = TRUE)
  signalHsmmNOHOM90 <- train_hsmm(seqs[seq90f], aaaggregation)
  
  list(signalHsmm = signalHsmm,
       signalHsmmNODEG = signalHsmmNODEG,
       signalHsmmNOHOM50 = signalHsmmNOHOM50,
       signalHsmmNOHOM90 = signalHsmmNOHOM90)
}
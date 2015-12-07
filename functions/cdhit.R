#' Filter sequrences using cd-hit
#'
#' Filters sequences using external software cd-hit.
#' @param input_seq list of input sequences
#' @param threshold threshold value
#' @param only_signal if \code{TRUE}, only signal peptides are filtered.
#' @param cdhit_path \code{character} path to cd-hit
#' @return vector of names of filtered sequences

cdhit <- function(input_seq, thresh = 0.5, word_length = 2, only_signal = TRUE, roi_length = 70,
                  cdhit_path = "/home/michal/cd-hit") {
  
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  cdhit <- paste0(cdhit_path, "/cd-hit -i ", input,  " -o ", output, " -c ", thresh, " -n ", word_length)
  
  if(only_signal) {
    #doesn't really make a lot of sense in the context of negative data set
    write.fasta(lapply(input_seq, function(single_seq) single_seq[1L:attr(single_seq, "signal")[2]]), 
                names = names(input_seq), input)
  } else {
    write.fasta(lapply(input_seq, function(single_seq) single_seq[1L:ifelse(length(single_seq) > roi_length, 
                                                                            roi_length, length(single_seq))]), 
                names = names(input_seq), input)
  }
  system(cdhit)
  res <- read.fasta(output)
  file.remove(input, output, paste0(output, ".clstr"))
  names(res)
}


#require(msa)
# require(signalHsmm)
# require(seqinr)
# plas_seqs <- read_uniprot("./plasmodium_benchmark_data/plas.txt", ft_names = "signal")
# 
# #write only SPs
# 
# 
# plas_fasta <- readAAStringSet("./plasmodium_benchmark_data/plas.fasta")
# plas_aln <- msaClustalW(plas_fasta, verbose = TRUE)

#' Filter sequrences using cd-hit
#'
#' Filters sequences using external software cd-hit.
#' @param input_seq list of input sequences
#' @param threshold threshold value
#' @param only_signal if \code{TRUE}, only signal peptides are filtered.
#' @return vector of names of filtered sequences

cdhit <- function(input_seq, thresh = 0.9, only_signal = TRUE) {
  
  input <- tempfile(tmpdir = getwd())
  output <- tempfile(tmpdir = getwd())
  cdhit <- paste0("/home/michal/cd-hit/cd-hit -i ", input,  " -o ", output, " -c ", thresh)
  
  if(only_signal) {
    write.fasta(lapply(input_seq, function(single_seq) single_seq[1L:attr(single_seq, "signal")[2]]), 
                names = names(input_seq), input)
  } else {
    write.fasta(input_seq, names = names(input_seq), input)
  }
  system(cdhit)
  res <- read.fasta(output)
  file.remove(input, output, paste0(output, ".clstr"))
  names(res)
}


require(signalHsmm)
require(seqinr)

plas_seqs <- read_uniprot("./plasmodium_benchmark_data/plas.txt", ft_names = "signal")
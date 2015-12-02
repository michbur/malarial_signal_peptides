sp_seqs <- read.fasta("./benchmark_data/sp2010_2015.fasta", seqtype = "AA")
nsp_seqs <- read.fasta("./benchmark_data/nsp2010_2015.fasta", seqtype = "AA")

sp_seqsf <- cdhit(sp_seqs, thresh = 0.5, word_length = 2, only_signal = FALSE)
nsp_seqsf <- cdhit(nsp_seqs, thresh = 0.5, word_length = 2, only_signal = FALSE)

set.seed(1)
chosen_negatives <- nsp_seqsf[sample(1L:length(nsp_seqsf), length(sp_seqsf), replace = FALSE)]
seq_NOHOM <- c(sp_seqs[sp_seqsf], nsp_seqs[chosen_negatives])

write.fasta(lapply(seq_NOHOM, function(i) i[1L:ifelse(length(i) > 150, 150, length(i))]), 
            names = names(seq_NOHOM), 
            file.out = "./benchmark_data/benchmark_data_NOHOM.fasta")

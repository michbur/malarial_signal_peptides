require(msa)
require(signalHsmm)
require(seqinr)

plas_seqs <- read_uniprot("./plasmodium_benchmark_data/plas.txt", ft_names = "signal")

#write only SPs
write.fasta(lapply(plas_seqs, function(single_seq) single_seq[1L:attr(single_seq, "signal")[2]]), 
            names = names(plas_seqs), "./plasmodium_benchmark_data/plas.fasta")

plas_fasta <- readAAStringSet("./plasmodium_benchmark_data/plas.fasta")
plas_aln <- msaClustalW(plas_fasta, verbose = TRUE)

#/home/michal/cd-hit/cd-hit -i plas.fasta -o plas_filtered.fasta


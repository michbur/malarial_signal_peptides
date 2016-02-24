


get_aa_freqs <- function(seqs, enc, taxon) {
  nondeg_sig <- do.call(rbind, lapply(seqs, function(single_seq) {
    aas <- single_seq[2L:attr(single_seq, "signal")[2]]
    as.matrix(count_ngrams(matrix(aas, nrow = 1), 1, a()[-1]))/length(aas)
  }))
  colnames(nondeg_sig) <- unlist(strsplit(colnames(nondeg_sig), "_0"))
  
  nondeg_mature <- do.call(rbind, lapply(seqs, function(single_seq) {
    aas <- single_seq[(attr(single_seq, "signal")[2] + 1):length(single_seq)]
    as.matrix(count_ngrams(matrix(aas, nrow = 1), 1, a()[-1]))/length(aas)
  }))
  colnames(nondeg_mature) <- unlist(strsplit(colnames(nondeg_mature), "_0"))
  
  deg_sig <- do.call(rbind, lapply(seqs, function(single_seq) {
    aas <- degenerate(tolower(single_seq[2L:attr(single_seq, "signal")[2]]),
                      enc)
    as.matrix(count_ngrams(matrix(aas, nrow = 1), 1, 1L:4))/length(aas)
  }))
  
  colnames(deg_sig) <- unlist(strsplit(colnames(deg_sig), "_0"))
  
  deg_mature <- do.call(rbind, lapply(seqs, function(single_seq) {
    aas <- degenerate(tolower(single_seq[(attr(single_seq, "signal")[2] + 1):length(single_seq)]), enc)
    as.matrix(count_ngrams(matrix(aas, nrow = 1), 1, 1L:4))/length(aas)
  }))
  colnames(deg_mature) <- unlist(strsplit(colnames(deg_mature), "_0"))
  
  list(deg = data.frame(taxon = taxon,
                        rbind(data.frame(type = "signal", deg_sig),
                              data.frame(type = "mature", deg_mature))),
       nondeg = data.frame(taxon = taxon,
                           rbind(data.frame(type = "signal", nondeg_sig),
                                 data.frame(type = "mature", nondeg_mature))))
}


freq_other <- get_aa_freqs(read_uniprot(paste0(pathway, "signal_peptides.txt"), ft_names = "signal"), 
                           enc_region[["best_sens_raw"]], "other")

freq_plas <- get_aa_freqs(read_uniprot("./plasmodium_benchmark_data/plas.txt", ft_names = "signal"), 
                           enc_region[["best_sens_raw"]], "plasmodium")

freq_deg <- rbind(freq_other[["deg"]], freq_plas[["deg"]])
write.csv2(freq_deg, file = "./frequency_analysis/freq_deg.csv")

freq_nondeg <- rbind(freq_other[["nondeg"]], freq_plas[["nondeg"]])
write.csv2(freq_nondeg, file = "./frequency_analysis/freq_nondeg.csv")

mfreq_deg <- melt(freq_deg)
mfreq_nondeg <- melt(freq_nondeg)

ggplot(mfreq_deg, aes(x = variable, y = value, fill = taxon)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(taxon ~ type)

ggplot(mfreq_nondeg, aes(x = variable, y = value, fill = taxon)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(taxon ~ type)

library(ica)
ica_deg <- icafast(freq_deg[, -c(1L:2)], nc = 2)

ica_deg_df <- cbind(freq_deg[, c(1L:2)], ica_deg[["S"]])

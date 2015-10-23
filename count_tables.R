library(signalHsmm)
library(seqinr)
library(reshape2)
library(dplyr)
library(biogram)

np_seq <- read_uniprot("non_plas.txt", ft_names = "signal")
p_seq <- read_uniprot("plas.txt", ft_names = "signal")

#signal: if TRUE, signal peptide, if FALSE, mature protein
get_freqs <- function(list_of_prots) {
  aa <- a()[-1]
  count_list <- lapply(list_of_prots, function(i) {
    end <- attr(i, "signal")[2]
    #signal peptide
    sp <- seqinr::count(i[1L:end], 1, alphabet = a()[-1])
    #signal peptide freq
    sp_freq <- sp/sum(sp)
    #mature protein
    mp <- seqinr::count(i[(end + 1):length(i)], 1, alphabet = a()[-1])
    mp_freq <- mp/sum(mp)
    sp_nfreq <- sp_freq/mp_freq
    data.frame(what = c(rep("mature protein", 2), rep("signal peptide", 3)),
               type = c(rep(c("count", "frequency"), 2), "normalized frequency"),
               rbind(mp, mp_freq, sp, sp_freq, sp_nfreq))
    
  })
  count_df <- do.call(rbind, count_list)
  count_df <- data.frame(prot_name = sapply(strsplit(rownames(count_df), ".", fixed = TRUE), first),
                         count_df)
  rownames(count_df) <- NULL
  count_df
}
  
whole_data <- rbind(data.frame(plasmodium = "no", get_freqs(np_seq)),
                    data.frame(plasmodium = "yes", get_freqs(p_seq)))
write.csv2(whole_data, file = "sp_freqs.csv", row.names = FALSE)
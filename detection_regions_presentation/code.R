library(ggplot2)
library(pbapply)
library(seqinr)
library(dplyr)
library(reshape2)
library(signalHsmm)
library(biogram)

if(Sys.info()["nodename"] == "phobos" )
  pathway <- "/home/michal/Dropbox/signal-peptide2_data/"

if(Sys.info()["nodename"] == "MICHALKOMP" )
  pathway <- "C:/Users/Michal/Dropbox/signal-peptide2_data/"


group2df <- function(group_list, caption = NULL, label = NULL) {
  data.frame(ID = 1L:length(group_list), 
             Groups = sapply(group_list, function(i)
               paste0(toupper(sort(i)), collapse = ", ")))
}

degreg <-  function(single_encoding, nhc_borders, pos_seqs) {
  single_encoding <- lapply(single_encoding, toupper)
  res <- lapply(1L:length(pos_seqs), function(seq_id) {
    single_seq <- pos_seqs[[seq_id]]
    aa_n <- single_seq[1L:(nhc_borders[seq_id, "start_h"] - 1)]
    aa_h <- single_seq[nhc_borders[seq_id, "start_h"]:(nhc_borders[seq_id, "start_c"] - 1)]
    aa_c <- single_seq[nhc_borders[seq_id, "start_c"]:(nhc_borders[seq_id, "cs"])]
    aa_mature <- single_seq[(nhc_borders[seq_id, "cs"] + 1):length(single_seq)]
    res <- list(aa_n = aa_n,
                aa_h = aa_h,
                aa_c = aa_c,
                aa_m = aa_mature)
    lapply(res, function(i) 
      degenerate(i, single_encoding))
  }) 
  tab <- data.frame(region = unlist(lapply(c("n", "h", "c", "m"), function(i) rep(i, length(single_encoding)))), 
                    do.call(rbind, lapply(1L:4, function(single_region) {
                      res <- factor(unlist(lapply(res, function(single_seq)
                        single_seq[[single_region]])), levels = as.character(1L:length(single_encoding)))
                      res_tab <- data.frame(table(res))
                      cbind(res_tab, prop = res_tab[["Freq"]]/length(res))
                    })))
  colnames(tab) <- c("region", "group", "count", "freq")
  tab[["region"]] <- factor(tab[["region"]], levels = c("n", "h", "c", "m"))
  levels(tab[["region"]]) <- c("n-region", "h-region", "c-region", "Mature\nprotein")
  levels(tab[["group"]]) <- paste0("Group ", levels(tab[["group"]]))
  
  tab
}


pos_seqs <- read_uniprot(paste0(pathway, "signal_peptides.txt"), ft_names = "signal")
longer_seqs <- pos_seqs[lengths(pos_seqs) > 100]
nhc_borders <- t(sapply(longer_seqs, find_nhc))

sumcleaves <- do.call(rbind, lapply(longer_seqs, function(single_seq) {
  sig <- attr(single_seq, "signal")[2]
  single_seq[(sig - 4):(sig + 4)]
})) %>% melt(value.name = "aa") %>% 
  rename(prot = Var1, pos = Var2) %>% 
  group_by(pos, aa) %>% 
  summarise(n = length(aa)) %>%
  ungroup %>%
  mutate(pos = pos - 6)

ggplot(sumcleaves, aes(x = aa, y = n, fill = aa)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ pos)


cs_aa <- c("A", "G", "L", "V", "P")

cs_enc <- list(`1` = cs_aa,
               `2` = setdiff(a()[-1], cs_aa))

cs_dat <- degreg(cs_enc, nhc_borders, longer_seqs)

ggplot(cs_dat, aes(x = group, y = freq, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ region, nrow = 1) + 
  scale_y_continuous("Frequency") +
  scale_x_discrete("Region\n")


nregions <- do.call(rbind, lapply(longer_seqs, function(single_seq) {
  single_seq[2:7]
})) %>% melt(value.name = "aa") %>% 
  rename(prot = Var1, pos = Var2) %>% 
  group_by(pos, aa) %>% 
  summarise(n = length(aa)) %>%
  ungroup 
  
ggplot(nregions, aes(x = aa, y = n, fill = aa)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ pos)

nreg_enc <- list(`1` = c("A", "G", "L", "V"), #both n-region and cs
                 `2` = c("G", "P", "K"), #only cs
                 `3` = c("R", "S", "T"), #only n-region
                 `4` = c("C", "D", "E", "F", "H", "I", "M", "N", "Q", "W", "Y"))

nreg_dat <- degreg(nreg_enc, nhc_borders, longer_seqs)
ggplot(nreg_dat, aes(x = group, y = freq, fill = group)) +
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ region, nrow = 1) + 
  scale_y_continuous("Frequency") +
  scale_x_discrete("Region\n")
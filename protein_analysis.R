preds <- read.csv("./plasmodium_protein_analysis/preds.csv") > 0.5
ets <- c(rep(TRUE, 51), rep(FALSE, 211))

cs <- c(18, 20, 25, 20, 16, 21, 24, 27, 23, 16, 23, 19, 24, 34, 23,
        16, 23, 20, 21, 22, 16, 24, 21, 19, 16, 21, 34, 27, 30, 20, 21,
        41, 21, 24, 26, 22, 29, 33, 20, 26, 21, 23, 22, 24, 21, 23, 22,
        25, 32, 29, 25)

all_TP <- data.frame(preds) %>% 
  slice(1L:51) %>%
  colSums %>%
  data.frame() %>% 
  select(TP = 1) %>%
  mutate(classifier = rownames(.),
         signalHsmm = grepl("signalHsmm", classifier))

data.frame(preds) %>% 
  slice(1L:51) %>%
  select(grep("_87", colnames(.)))

only_sig <- data.frame(preds) %>% 
  slice(1L:51) %>%
  select(signalHsmmNOHOM50_10, signalP3nn) %>% 
  rowSums == 1 

all_prots <- read.fasta("./plasmodium_benchmark_data/benchmark_plas_data_NOHOM.fasta",
                        seqtype = "AA")[1L:51]

get_signals <- function(all_prots, cs, u) {
  do.call(rbind, lapply(1L:length(all_prots), function(id) {
    table(factor(degenerate(all_prots[[id]][2L:cs[id]], u), levels = names(u)))/(cs[id] - 1)
    })) %>% 
    data.frame %>% 
    mutate(only_sig = only_sig, prot_id = paste0("Protein: ", 1L:nrow(.))) %>% 
    melt(variable.name = "aa")
}

aa <- a()[-1]
names(aa) <- a()[-1]
signals <- get_signals(all_prots, cs, aa)
signals_deg <- get_signals(all_prots, cs, sapply(enc_region[["best_sens_raw"]], toupper))

ggplot(filter(signals, only_sig), aes(x = aa, y = value)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ prot_id) +
  scale_x_discrete("Amino acid") +
  scale_y_continuous("Frequency") +
  my_theme

ggplot(filter(signals_deg, only_sig), aes(x = aa, y = value)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ prot_id) +
  scale_x_discrete("Amino acid") +
  scale_y_continuous("Frequency") +
  my_theme


deg_freq_plot <- function(x) {
  agg_signals <- group_by(x, only_sig, aa) %>% 
    summarise(value = mean(value)) %>% 
    ungroup
  
  ggplot(agg_signals, aes(x = aa, y = value, fill = !only_sig)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete("Amino acid") +
    scale_y_continuous("Mean frequency") +
    scale_fill_discrete("Recognized by signalP") +
    my_theme
}

cs_dat <- data.frame(cs, only_sig)

save(all_TP, signals, signals_deg, cs_dat, my_theme, file = "./plasmodium_protein_analysis/presentation.RData")

only_sig <- data.frame(preds) %>% 
  slice(1L:51) %>%
  select(signalHsmmNOHOM50_10, signalP3nn) %>% 
  rowSums == 1 
only_sig_sp4 <- data.frame(preds) %>% 
  slice(1L:51) %>%
  select(signalHsmmNOHOM50_10, signalP41notm) %>% 
  rowSums == 1 
dput(names(all_prots)[only_sig])
names(all_prots)[only_sig_sp4] 

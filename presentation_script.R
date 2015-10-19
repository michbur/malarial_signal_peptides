library(signalHsmm)
library(seqinr)
library(reshape2)
library(dplyr)
library(biogram)

#taxonomy:"Eukaryota [2759]" NOT taxonomy:"Plasmodium [5820]" annotation:(type:signal evidence:ECO_0000269) AND reviewed:yes
np_seq <- read_uniprot("non_plas.txt", ft_names = "signal")

#taxonomy:"Plasmodium [5820]" annotation:(type:signal evidence:manual) AND reviewed:yes
p_seq <- read_uniprot("plas.txt", ft_names = "signal")



get_signals <- function(list_of_prots) {
  aa <- a()[-1]
  lapply(list_of_prots, function(i) {
    end <- attr(i, "signal")[2]
    sig_freq <- data.frame(sig = TRUE, table(factor(i[1L:end], level = aa)))
    nsig_freq <- data.frame(sig = FALSE, table(factor(i[(end + 1):length(i)], level = aa)))
    nsig_freq[["Freq"]] <- nsig_freq[["Freq"]]/sum(nsig_freq[["Freq"]])
    res <- rbind(sig_freq, nsig_freq)
    total_freq <- data.frame(table(factor(i, levels = aa)))[["Freq"]]
    #res[["Freq"]] <- abs((res[["Freq"]] - total_freq))
    res
  })
}

np_freq <- do.call(rbind, get_signals(np_seq))
np_freq[["prot"]] <- sapply(strsplit(rownames(np_freq), ".", fixed = TRUE), first)
np_freq[["plasmodium"]] <- FALSE

p_freq <- do.call(rbind, get_signals(p_seq))
p_freq[["prot"]] <- sapply(strsplit(rownames(p_freq), ".", fixed = TRUE), first)
p_freq[["plasmodium"]] <- TRUE

freq_dat <- rbind(np_freq, p_freq)
colnames(freq_dat) <- c("sig", "aa", "freq", "prot", "plasmodium")


quick_diff <- function(x)
  abs(Reduce("-", x))
diff_dat <- group_by(freq_dat, sig, aa, plasmodium) %>% 
  summarise(mfreq = mean(freq)) %>%
  summarise(dmfreq = quick_diff(mfreq)) %>%
  filter(dmfreq > quantile(dmfreq, 0.75)) 

mean_diff_dat <- group_by(freq_dat, sig, aa, plasmodium) %>% 
  summarise(mfreq = mean(freq)) %>%
  filter(aa %in% diff_dat[["aa"]])


get_degs<- function(list_of_prots) {
  aa <- as.character(1L:4)
  lapply(list_of_prots, function(i) {
    end <- attr(i, "signal")[2]
    sig_freq <- data.frame(sig = TRUE, 
                           table(factor(degenerate(i[1L:end], aaaggregation), level = aa)))
    nsig_freq <- data.frame(sig = FALSE, 
                            table(factor(degenerate(i[(end + 1):length(i)], aaaggregation), level = aa)))
    nsig_freq[["Freq"]] <- nsig_freq[["Freq"]]/sum(nsig_freq[["Freq"]])
    res <- rbind(sig_freq, nsig_freq)
    total_freq <- data.frame(table(factor(i, levels = aa)))[["Freq"]]
    #res[["Freq"]] <- abs((res[["Freq"]] - total_freq))
    res
  })
}

np_freq <- do.call(rbind, get_degs(np_seq))
np_freq[["prot"]] <- sapply(strsplit(rownames(np_freq), ".", fixed = TRUE), first)
np_freq[["plasmodium"]] <- FALSE

p_freq <- do.call(rbind, get_degs(p_seq))
p_freq[["prot"]] <- sapply(strsplit(rownames(p_freq), ".", fixed = TRUE), first)
p_freq[["plasmodium"]] <- TRUE

freq_dat <- rbind(np_freq, p_freq)
colnames(freq_dat) <- c("sig", "aa", "freq", "prot", "plasmodium")

mean_group_dat <- group_by(freq_dat, sig, aa, plasmodium) %>% 
  summarise(mfreq = mean(freq))

length_dat <- melt(rbind(data.frame(plasmodium = FALSE, 
                               sig = sapply(np_seq, function(i) attr(i, "sig")[2]),
                               len = lengths(np_seq)),
                    data.frame(plasmodium = TRUE, 
                               sig = sapply(p_seq, function(i) attr(i, "sig")[2]),
                               len = lengths(p_seq))))

filter(length_dat, variable == "len") %>% 
  group_by(plasmodium) %>%
  summarize(mean = mean(value), median = median(value))


get_regions_freq <- function(list_of_prots) {
  aa <- as.character(1L:4)
  lapply(list_of_prots, function(i) {
    borders <- find_nhc(i, attr(i, "signal"))
    nr <- data.frame(region = "n", table(factor(degenerate(i[1L:borders[2]], aaaggregation))))
    hr <- data.frame(region = "h", table(factor(degenerate(i[(borders[2]+1):borders[3]], aaaggregation))))
    cr <- data.frame(region = "c", table(factor(degenerate(i[(borders[3]+1):borders[4]], aaaggregation))))
    rbind(nr, hr, cr)
  })
}


np_freq <- do.call(rbind, get_regions_freq(np_seq))
np_freq[["prot"]] <- sapply(strsplit(rownames(np_freq), ".", fixed = TRUE), first)
np_freq[["plasmodium"]] <- FALSE

p_freq <- do.call(rbind, get_regions_freq(p_seq))
p_freq[["prot"]] <- sapply(strsplit(rownames(p_freq), ".", fixed = TRUE), first)
p_freq[["plasmodium"]] <- TRUE

freq_dat <- rbind(np_freq, p_freq)
colnames(freq_dat) <- c("region", "aa", "freq", "prot", "plasmodium")

region_dat <- group_by(freq_dat, region, aa, plasmodium) %>% 
  summarise(mfreq = mean(freq))

get_region_length <- function(list_of_prots) {
  dat <- data.frame(t(sapply(list_of_prots, function(i) find_nhc(i, attr(i, "signal")))))
  mutate(dat, len_n = start_h - start_n,
         len_h = start_c - start_h,
         len_c = cs - start_c) %>%
    select(len_n, len_h, len_c) %>%
    melt
}

reg_length <- rbind(data.frame(plasmodium = FALSE, get_region_length(np_seq)),
                    data.frame(plasmodium = TRUE, get_region_length(p_seq)))



ggplot(filter(reg_length, !plasmodium), aes(x = value, fill = plasmodium)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ variable)

save(diff_dat, mean_diff_dat, mean_group_dat, length_dat, region_dat,
     reg_length, file = "apicomplexa_signalpeptides.RData")

seq50_10 <- read_uniprot("./training_data/sp1950_2010.txt", ft_names = "signal")

t(sapply(seq50_10, find_nhc))

find_nhc4 <- function(protein, signal = NULL) {
    protein <- toupper(protein)
    if (is.null(signal)) 
      signal <- attr(protein, "signal")
    
    sig <- protein[signal[1]:signal[2]]
    
    start_c <- signal[2] - 4
    
    start_n <- 1
    start_h <- 2
    noh <- 0
    
    while(noh < 3 && start_h < start_c) {
      start_h <- start_h + 1
      noh <- noh + ifelse(sig[start_h] %in% c("F", "I", "L", "V", "W"), 1, 0)
      noh <- ifelse(noh < 1, 0, noh)
    }
    
    start_h <- start_h - 3
    
    c(start_n = 1, start_h = start_h, start_c = start_c, cs = signal[2])
}




find_nhc3 <- function(protein, signal = NULL) {
  protein <- toupper(protein)
  if (is.null(signal)) 
    signal <- attr(protein, "signal")
  
  protein[(signal[1] + 1):(signal[1] + 9)]
}

sigs <- sapply(seq50_10, function(protein)
  attr(protein, "signal")[2])

sapply(seq50_10[sigs == 19], find_nhc3) %>% t %>% 
  apply(2, function(i) table(factor(i, levels = a()[-1]))) %>%
  melt %>%
  ggplot(aes(x = Var1, y = value, fill = Var1)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Var2)

t(sapply(seq50_10, find_nhc4))

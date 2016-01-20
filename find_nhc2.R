seq50_10 <- read_uniprot("./training_data/sp1950_2010.txt", ft_names = "signal")

t(sapply(seq50_10, find_nhc))

find_nhc2 <- function(protein, signal = NULL) {
    protein <- toupper(protein)
    if (is.null(signal)) 
      signal <- attr(protein, "signal")
    
    sig <- protein[signal[1]:signal[2]]
    
    start_c <- signal[2] - 4
    
    start_h <- start_c - 6
    #if statement to prevent negative values
    if (start_h > 1) {
      #nonh number of nonhydrophobic residues
      nonh <- 0
      #noc number of charged
      noc <- 0
      while(nonh < 3 && noc == 0 && start_h > 1) {
        start_h <- start_h - 1
        nonh <- nonh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), -1, 1)
        nonh <- ifelse(nonh < 0, 0, nonh)
        noc <- ifelse(sig[start_h] %in% c("R", "H", "K", "D", "E"), 1, 0)
      }
    } else {
      start_h <- 1
    }
    
    prestart_c <- start_c - 1
    noh <- 0
    while(noh == 0 && start_h < prestart_c) {
      start_h <- start_h + 1
      noh <- noh + ifelse(sig[start_h] %in% c("A", "I", "L", "M", "F", "W", "V"), 1, 0)
    }
    #c(start_n = signal[1], start_h = start_h, start_c = start_c, cs = signal[2])
    c(start_n = 1, start_h = start_h, start_c = start_c, cs = signal[2])
}

t(sapply(seq50_10, find_nhc2))

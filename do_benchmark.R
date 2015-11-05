

# library(ggplot2)
# dat <- data.frame(real = real_labels, pred = signalHsmm_preds)
# 
# ggplot(dat, aes(x = as.factor(real), y = pred)) +
#   geom_boxplot()
# 
# 
# library(dplyr)
# group_by(dat, real) %>% summarise(median(pred))

# source of proteins in signalHsmm1986
# prot86 <- read_uniprot("./training_data/sp1950_1987.txt", ft_names = "signal")
# table(unname(sapply(prot86, function(i) attr(i, "OS"))))
dat <- do.call(rbind, strsplit(readLines("./plasmodium_benchmark/signalP30.txt")[-1], " +"))[, c(3, 13, 17, 19)]
res <- t(apply(dat, 1, as.numeric))
#predictions for signalP 3.0
preds <- round(res[, c(2, 4)], 2) > 0.5
colnames(preds) <- c("signalP3nn", "signalP3hmm")

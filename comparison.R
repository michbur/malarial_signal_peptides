read_signalp41 <- function(connection) {
  all_lines <- readLines(connection)[-c(1L:2)]
  do.call(rbind, lapply(all_lines, function(i) {
    line <- strsplit(i, " ")[[1]]
    line <- line[!line == ""]
    res <- data.frame(sp.probability = line[10] == "Y",
                      sp.start = ifelse(line[10] == "Y", 1, NA),
                      sp.end = ifelse(line[10] == "Y", as.numeric(line[5]) - 1, NA))
    rownames(res) <- line[1]
    res
  }))
}

read_predsi <- function(connection) {
  dat <- read.table(connection, sep = "\t")
  data.frame(sp.probability = dat[, 4] == "Y",
             sig.start = ifelse(dat[, 4] == "Y", 1, NA),
             sig.end = ifelse(dat[, 4] == "Y", as.numeric(dat[, 3]), NA),
             row.names = dat[, 1])
}

read_phobius <- function(connection) {
  all_lines <- readLines(connection)
  all_lines <- all_lines[-1]
  splited <- strsplit(all_lines, " ")
  #remove "" characters
  purged <- t(sapply(splited, function(i) i[i != ""]))
  cl_sites <- sapply(which(purged[, 3] == "Y"), function(i)
    as.numeric(strsplit(strsplit(purged[i,4], "/")[[1]][1], "c")[[1]][[2]]))
  res <- data.frame(sp.probability = purged[, 3] == "Y",
                    sig.start = ifelse(purged[, 3] == "Y", 1, NA),
                    sig.end = rep(NA, nrow(purged)), 
                    row.names = purged[, 1])
  res[purged[, 3] == "Y", "sig.end"] <- cl_sites
  res
}

read_philius <- function(connection) {
  require(XML)
  all_dat <- xmlToList(xmlTreeParse(connection, asTree = TRUE))
  seq_dat_id <- 1L:(length(all_dat)/2)*2
  #data for table
  table_dat <- sapply(seq_dat_id, function(i) 
    unlist(all_dat[i][[1]][[1]][c(24, 22)]))
  cleaved <- sapply(table_dat, function(i)
    !(is.null(i[1]) || is.na(i[1])))
  res <- data.frame(sp.probability = cleaved,
                    sig.start = ifelse(cleaved, 1, NA),
                    sig.end = rep(NA, length(seq_dat_id)),
                    row.names = unlist(all_dat[1L:(length(all_dat)/2)*2 - 1]))
  res[cleaved, "sig.end"] <- as.numeric(sapply(table_dat[cleaved], function(i) i[2])) - 1
  res
}

library(signalHsmm)
library(seqinr)
library(hmeasure)


signal.hsmm1987_preds <- pred2df(predict(signal.hsmm1987, c(plas_pos, plas_neg)))
signal.hsmm_preds <- pred2df(run_signalHsmm(c(plas_pos, plas_neg)))

real_labels <- c(rep(1, length(plas_pos)), rep(0, length(plas_neg)))

metrics <- do.call(rbind, lapply(list(signalPnotm = read_signalp41("./benchmark/signaP41notm.txt"), 
                                      signalPtm = read_signalp41("./benchmark/signaP41tm.txt"), 
                                      predsi = read_predsi("./benchmark/predsi.txt"),
                                      phobius = read_phobius("./benchmark/phobius.txt"),
                                      philius = read_philius("./benchmark/philius.xml"),
                                      signalHsmm = signal.hsmm1987_preds,
                                      signalHsmm1987 = signal.hsmm1987_preds), function(predictor)
                                        HMeasure(real_labels, predictor[["sp.probability"]])[["metrics"]]))
TP <- as.numeric(metrics[["TP"]])
FP <- as.numeric(metrics[["FP"]])
TN <- as.numeric(metrics[["TN"]])
FN <- as.numeric(metrics[["FN"]])


write.table(round(cbind(metrics[, c("AUC", "H")], 
                  MCC = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))), 4), sep = "\t", file = "grant_short.txt")


other_software <- t(do.call(rbind, lapply(list(signalPnotm = read_signalp41("./benchmark/signaP41notm.txt"), 
                                      signalPtm = read_signalp41("./benchmark/signaP41tm.txt"), 
                                      predsi = read_predsi("./benchmark/predsi.txt"),
                                      phobius = read_phobius("./benchmark/phobius.txt")), function(predictor)
                                        predictor[["sp.probability"]])))
write.csv2(other_software, file = "other_soft.csv")

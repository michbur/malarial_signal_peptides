#taxonomy:"Plasmodium [5820]" annotation:(type:signal evidence:manual) AND reviewed:yes
#111 proteins
#taxonomy:"Plasmodium [5820]" NOT annotation:(type:signal evidence:manual) AND reviewed:yes
#361 proteins


# library(seqinr)
# dat <- read.fasta("benchmark_plas_data.fasta")
# #truncated data for signalP 3.0
# short_dat <- lapply(dat, function(i) i[1L:200])
# write.fasta(short_dat, names = names(dat), file.out = "benchmark_plas_data_tr.fasta")

# library(hmeasure)
# real_labels <- c(rep(1, 102), rep(0, 358))
# 
# HMeasure(real_labels, read.csv2("./plasmodium_benchmark/benchmark_plas_other.csv"))[["metrics"]]


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

read_signalp3 <- function(connection) {
  #columns: nn Ymax, nn Ymax position, hmm Ymax 
  dat <- do.call(rbind, strsplit(readLines(connection)[-1], " +"))[, c(5, 6, 17, 19)]
  res <- data.frame(t(apply(dat, 1, as.numeric)))
  colnames(res) <- c("nn_prob", "nn_pos", "hmm_prob", "hmm_pos")
  res
}

read_signalp3nn <- function(connection) {
  res <- read_signalp3(connection)
  data.frame(sp.probability = res[["nn_prob"]] > 0.5, sp.start = 1, sp.end = res[["nn_pos"]])
}

read_signalp3hmm <- function(connection) {
  res <- read_signalp3(connection)
  data.frame(sp.probability = res[["hmm_prob"]] > 0.5, sp.start = 1, sp.end = res[["hmm_pos"]])
}

read_signalp3("plasmodium_benchmark_results/signalP30.txt")

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

signal.hsmm1987_preds <- pred2df(predict(signal.hsmm1987, c(plas_pos, plas_neg)))
signal.hsmm_preds <- pred2df(run_signalHsmm(c(plas_pos, plas_neg)))

real_labels <- c(rep(1, length(plas_pos)), rep(0, length(plas_neg)))

metrics <- do.call(rbind, lapply(list(signalP41notm = read_signalp41("./plasmodium_benchmark_results/signaP41notm.txt"), 
                                      signalP41tm = read_signalp41("./plasmodium_benchmark_results/signaP41tm.txt"), 
                                      signalP3nn = read_signalp3nn("plasmodium_benchmark_results/signalP30.txt"),
                                      signalP3hmm = read_signalp3hmm("plasmodium_benchmark_results/signalP30.txt"),
                                      predsi = read_predsi("./plasmodium_benchmark_results/predsi.txt"),
                                      phobius = read_phobius("./plasmodium_benchmark_results/phobius.txt"),
                                      philius = read_philius("./plasmodium_benchmark_results/philius.xml"),
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
#write.csv2(other_software, file = "other_soft.csv")

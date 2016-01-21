#' Read prediction of other software
#'
#' Read predictions from other signal peptide predictors provided that they are in 
#' the proper format and in specified directory.
#'
#' @param directory_name for example "./plasmodium_benchmark_results"
#'
#' @return a data frame of predictions

read_other_software <- function(directory_name) {
  
  read_signalp41 <- function(connection) {
    dat <- read.table(connection)
    
    data.frame(#sp.probability = dat[, 10] == "Y",
      sp.probability = dat[, 9],
      sp.start = ifelse(dat[, 10] == "Y", 1, NA),
      sp.end = ifelse(dat[, 10] == "Y", as.numeric(dat[, 5]) - 1, NA))
  }
  
  read_signalp3 <- function(connection) {
    #columns: nn Ymax, nn Ymax position, hmm Ymax 
    dat <- do.call(rbind, strsplit(readLines(connection)[-1], " +"))[, c(13, 6, 17, 19)] #mean D for NN and Signal peptide probability for HMM
    res <- data.frame(t(apply(dat, 1, as.numeric)))
    colnames(res) <- c("nn_prob", "nn_pos", "hmm_pos", "hmm_prob")
    res
  }
  
  read_signalp3nn <- function(connection) {
    res <- read_signalp3(connection)
    #data.frame(sp.probability = res[["nn_prob"]] > 0.5, sp.start = 1, sp.end = res[["nn_pos"]])
    data.frame(sp.probability = res[["nn_prob"]], sp.start = 1, sp.end = res[["nn_pos"]])
  }
  
  read_signalp3hmm <- function(connection) {
    res <- read_signalp3(connection)
    #data.frame(sp.probability = res[["hmm_prob"]] > 0.5, sp.start = 1, sp.end = res[["hmm_pos"]])
    data.frame(sp.probability = res[["hmm_prob"]], sp.start = 1, sp.end = res[["hmm_pos"]])
  }
  
  read_predsi <- function(connection) {
    dat <- read.table(connection, sep = "\t")
    data.frame(sp.probability = dat[, 2],
               #sp.probability = dat[, 4] == "Y",
               sig.start = ifelse(dat[, 4] == "Y", 1, NA),
               sp.end = ifelse(dat[, 4] == "Y", as.numeric(dat[, 3]), NA),
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
                      sp.end = rep(NA, nrow(purged)), 
                      row.names = purged[, 1])
    res[purged[, 3] == "Y", "sp.end"] <- cl_sites
    res
  }
  
  read_philius <- function(connection) {
    all_dat <- xmlToList(xmlTreeParse(connection, asTree = TRUE))
    seq_dat_id <- 1L:(length(all_dat)/2)*2
    #data for table
    table_dat <- sapply(seq_dat_id, function(i) 
      unlist(all_dat[i][[1]][[1]][c(24, 22)]))
    cleaved <- sapply(table_dat, function(i)
      !(is.null(i[1]) || is.na(i[1])))
    res <- data.frame(sp.probability = cleaved,
                      sig.start = ifelse(cleaved, 1, NA),
                      sp.end = rep(NA, length(seq_dat_id)),
                      row.names = unlist(all_dat[1L:(length(all_dat)/2)*2 - 1]))
    res[cleaved, "sp.end"] <- as.numeric(sapply(table_dat[cleaved], function(i) i[2])) - 1
    res
  }
  

  prob <- do.call(cbind, lapply(list(signalP41notm = read_signalp41(paste0(directory_name, "/signaP41notm.txt")), 
                             signalP41tm = read_signalp41(paste0(directory_name, "/signaP41tm.txt")), 
                             signalP3nn = read_signalp3nn(paste0(directory_name, "/signalP30.txt")),
                             signalP3hmm = read_signalp3hmm(paste0(directory_name, "/signalP30.txt")),
                             predsi = read_predsi(paste0(directory_name, "/predsi.txt")),
                             philius = read_philius(paste0(directory_name, "/philius.xml")),
                             phobius = read_phobius(paste0(directory_name, "/phobius.txt"))), function(predictor)
                               predictor[["sp.probability"]]))
  pos <- do.call(cbind, lapply(list(signalP41notm = read_signalp41(paste0(directory_name, "/signaP41notm.txt")), 
                                     signalP41tm = read_signalp41(paste0(directory_name, "/signaP41tm.txt")), 
                                     signalP3nn = read_signalp3nn(paste0(directory_name, "/signalP30.txt")),
                                     signalP3hmm = read_signalp3hmm(paste0(directory_name, "/signalP30.txt")),
                                     predsi = read_predsi(paste0(directory_name, "/predsi.txt")),
                                     philius = read_philius(paste0(directory_name, "/philius.xml")),
                                     phobius = read_phobius(paste0(directory_name, "/phobius.txt"))), function(predictor)
                                       predictor[["sp.end"]]))
  
  #remove cs position information for improbable signal peptides
  pos[prob < 0.5] <- NA
  
  list(prob = prob, pos = pos)
}

#' Calculate performance metrics
#'
#' Calculate performance metrics as \code{\link[hmeasure]{HMeasure}} plus 
#' MCC.
#'
#' @param real_labels real labels of data.
#' @param preds data.frame of predictions, preferably named.
#' @param threshold cut-off.
#'
#' @return a data frame of performance measures

calc_metrics <- function(real_labels, preds, threshold = 0.5) {
  metrics <- HMeasure(true.class = real_labels, scores = preds, threshold = threshold)[["metrics"]]
  
  TP <- as.numeric(metrics[["TP"]])
  FP <- as.numeric(metrics[["FP"]])
  TN <- as.numeric(metrics[["TN"]])
  FN <- as.numeric(metrics[["FN"]])
  
  MCC = (TP*TN - FP*FN)/sqrt((TP + FP)*(TP + FN)*(TN + FP)*(TN + FN))
  data.frame(metrics, MCC)
}

get_signalHsmm_preds <- function(list_of_predictors, connection) {
  seqs <- read.fasta(connection, seqtype = "AA")
  prob <- sapply(list_of_predictors, function(single_predictor)
    pred2df(predict(single_predictor, seqs))[, "sp.probability"])
  pos <- sapply(list_of_predictors, function(single_predictor)
    pred2df(predict(single_predictor, seqs))[, "sp.end"])
  colnames(prob) <- names(list_of_predictors)
  colnames(pos) <- names(list_of_predictors)
  #   colnames(res) <- substitute(list_of_predictors) %>%
  #     deparse %>%
  #     paste0(collapse = "") %>% #paste too long calls
  #     strsplit(., "(", fixed = TRUE) %>% 
  #     unlist %>%
  #     nth(2) %>%
  #     strsplit(x = ., split = ")", fixed = TRUE) %>%
  #     unlist %>%
  #     strsplit(., ",[ ]+") %>%
  #     unlist
  
  #remove cs position information for improbable signal peptides
  pos[prob < 0.5] <- NA
  
  list(prob = prob, pos = pos)
}

#' Format benchmark table
#'
#' Format table returned by \code{calc_metrics}. Automatically adds citations, bolds 
#' highest measures and so on.
#'
#' @param x \code{data.frame} of benchmark data - output of \code{calc_metrics}.
#' @param caption caption.
#' @param label label.
#'
#' @return a \code{xtable} object.
format_bench_table <- function(x, caption, label) {
  bold_max <- function(x) {
    nx <- as.numeric(x)
    x[nx == max(nx)] <- paste0("\\textbf{", x[nx == max(nx)], "}")
    x
  }
  
  tab <- data.frame(x)[, c("AUC", "Sens", "Spec", "MCC")] %>%
    format(digits = 4) %>% mutate(AUC = bold_max(AUC),
                                  Sens = bold_max(Sens),
                                  Spec = bold_max(Spec),
                                  MCC = bold_max(MCC))
  
  colnames(tab) <- c("AUC", "Sensitivity", "Specificity", "MCC")
  rownames(tab) <- c("signalP 4.1 (no tm) \\cite{2011petersensignalp}",
                     "signalP 4.1 (tm) \\cite{2011petersensignalp}",
                     "signalP 3.0 (NN) \\cite{2004bendtsenimproved}",
                     "signalP 3.0 (HMM) \\cite{2004bendtsenimproved}",
                     "PrediSi \\cite{2004hillerpredisi}",
                     "Phobius \\cite{2004klla}",
                     "Philius \\cite{2008reynoldstransmembrane}",
                     #"signalHsmm-2010 (k-mer)", 
                     "signalHsmm-2010", "signalHsmm-2010 (raw aa)", "signalHsmm2010 (hom. 90\\%)", "signalHsmm2010 (hom. 50\\%)", 
                     "signalHsmm-1987", "signalHsmm-1987 (raw aa)", "signalHsmm1987 (hom. 90\\%)", "signalHsmm1987 (hom. 50\\%)")
  
  rws <- seq(1, nrow(tab) - 1, by = 2)
  col <- rep("\\rowcolor[gray]{0.85}", length(rws))
  
  res <- print(xtable(tab, caption = caption, label = label), include.rownames = TRUE, booktabs = TRUE,
               add.to.row = list(pos = as.list(rws), command = col), print.results = FALSE, caption.placement = "top",
               sanitize.text.function = identity, sanitize.rownames.function = identity)
  res
}

# count_signals <- function(connection) {
#   lines <- readLines(connection)
#   length(grep(pattern = "^FT   SIGNAL", lines))
# }

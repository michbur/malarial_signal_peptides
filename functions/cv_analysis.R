#' Get performance measures of cross-validation
#'
#' Get performance measures of cross-validation
#' @param folds cross-validation results (very large list)
#' @inheritParams calc_metrics
#'
#' @return a data frame of cross-validation performance measures

perf_rep <- function(folds, threshold = 0.5) {
  do.call(rbind, pblapply(1L:length(folds), function(repetition_id) {
    res <- t(sapply(folds[[repetition_id]], function(single_fold)
      rowMeans(sapply(single_fold, function(single_group) 
        #calc_metrics is defined in benchmark_functions.R
        unlist(calc_metrics(as.numeric(!is.na(single_group[, "cs_real"])), single_group[, "prob"], threshold))
      ))
    ))
    
    res <- melt(res)
    
    colnames(res) <- c("encoding", "measure", "value")
    
    res[["encoding"]] <- factor(res[["encoding"]])
    cbind(repetition = factor(rep(repetition_id, nrow(res))), res)
  }))
}

#' Create data for cvplot
#'
#' Create data for cvplot (cross-validation results for all encodings).
#' @param rep_res results cross-validation  as produced by \code{\link{perf_rep}}
#'
#' @return a data frame ready for p1 plot.

#compute results of cv analysis
create_cvplotdat <- function(rep_res) {
  mean_res <- rep_res %>% filter(measure %in% c("AUC", "Sens", "Spec", "MCC")) %>%
    group_by(encoding, measure) %>% summarise(mean_value = mean(value, na.rm = TRUE)) %>% ungroup
  
  best_enc <- lapply(c("AUC", "Sens", "Spec", "MCC"), function(i)
    mean_res %>% filter(measure == i) %>%
      filter(mean_value > quantile(mean_value, 0.9, na.rm = TRUE)) %>%
      droplevels %>%
      #arrange(desc(mean_value)) %>%
      select(encoding) %>%
      unlist %>%
      as.character %>%
      as.numeric)
  
  
  p1_dat <- rep_res %>% filter(#encoding %in% unique(unlist(best_enc)), 
    measure %in% c("AUC", "Spec", "Sens", "MCC")) %>% droplevels %>%
    group_by(encoding, measure) %>% summarize(mean = mean(value)) %>%
    dcast(encoding ~ measure, value.var = "mean")
  
  p1_dat <- p1_dat[!duplicated(p1_dat[, -1]), ]
  p1_dat[, "encoding"] <- rep("", nrow(p1_dat))
  p1_dat[p1_dat[, "Spec"] == max(p1_dat[, "Spec"]), "encoding"] <- "2"
  #p1_dat[p1_dat[, "Sens"] > 0.85 & p1_dat[, "Spec"] > 0.94, "encoding"] <- "3"
  p1_dat[p1_dat[, "Sens"] == max(p1_dat[, "Sens"]), "encoding"] <- "1"
  p1_dat[, "encoding"] <- as.factor(p1_dat[, "encoding"])
  p1_dat
}

#' Create p1 plot
#'
#' Create cross-validation plot and dynamically generate caption as well as table.
#' @param p1_dat processed data for the plot as produced by \code{\link{create_cvplotdat}}.
#'
#' @return a list, the first element is a plot, the second is a caption, the third is 
#' the table.

plot_cvplot <- function(p1_dat) {
  p1 <- ggplot(p1_dat, aes(x = Sens, y = Spec, label = encoding, colour = encoding == "", fill = encoding == "")) +
  geom_point(size = 5, shape = 21) +
  geom_text(size = 6, hjust = -0.75, vjust = 1) +
  scale_colour_manual(values = c("red","blue")) + 
  scale_fill_manual(values = c(adjustcolor("red", 0.25), adjustcolor("blue", 0.25))) + 
  scale_x_continuous("Sensitivity\n") +
  scale_y_continuous("Specificity") + 
  my_theme +
  guides(colour = FALSE, fill = FALSE)
  
  #caption for cvres
  cpt <- paste0("Results of cross-validation. 1. The encoding providing the best sensitivity (AUC = ", 
           round(mean(p1_dat[p1_dat[, "encoding"] == "1", "AUC"]), 4),
           ", MCC = ", 
           round(mean(p1_dat[p1_dat[, "encoding"] == "1", "MCC"]), 4),
           "). 2. The encoding providing the best specificity (AUC = ", 
           round(mean(p1_dat[p1_dat[, "encoding"] == "2", "AUC"]), 4),
           ", MCC = ", 
           round(mean(p1_dat[p1_dat[, "encoding"] == "2", "MCC"]), 4),
           ").")
  
  #cv table for the best encoding (69)
  tab <- filter(rep_res, encoding == 69, measure %in% c("AUC", "Sens", "Spec", "MCC")) %>%
    group_by(measure) %>% 
    summarize(m_value = mean(value), sd_value = sd(value)) %>% 
    droplevels %>% 
    data.frame
  
  levels(tab[["measure"]]) <- c("AUC", "Sensitivity", "Specificity", "MCC")
  
  rws <- seq(1, nrow(tab) - 1, by = 2)
  colnames(tab) <- c("Measure", "Mean", "SD")
  col <- rep("\\rowcolor[gray]{0.85}", length(rws))
  xtab <- print(xtable(tab, caption = "Performance measures for the best encoding. 60 repetitions of cross-validation.", 
               label = "tab:perfmeas", digits = 4), include.rownames = FALSE, booktabs = TRUE,
        add.to.row = list(pos = as.list(rws), command = col), print.results = FALSE,
        caption.placement = "top")
  
  
  list(plot = p1, 
       cpt = cpt, 
       xtab = xtab)
}

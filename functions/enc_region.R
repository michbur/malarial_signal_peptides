#' Create encoding-related plots and tables
#'
#' Create plots and tables related to the best performing encodings.
#' @param p1_dat cross-validation of encodings as returned by 
#' \code{\link{create_cvplotdat}}
#'
#' @return a list of four: a) plot of properties of amino acids in the groups of
#' encodings, b) frequency of groups in regions of signal peptide, c) the best 
#' sensitivity encoding and d) the best specificity encoding.

create_enc_region <- function(p1_dat) {
  
  group2df <- function(group_list, caption = NULL, label = NULL) {
    tab <- data.frame(Groups = sapply(group_list, function(i)
      paste0(toupper(sort(i)), collapse = ", ")))
    rws <- seq(1, nrow(tab) - 1, by = 2)
    col <- rep("\\rowcolor[gray]{0.85}", length(rws))
    print(xtable(tab, caption = caption, label = label), include.rownames = FALSE, booktabs = TRUE,
          add.to.row = list(pos = as.list(rws), command = col), print.results = FALSE, caption.placement = "top")
  }
  
  data(aaindex)
  
  #normalized values of amino acids 
  aa_nvals <- t(sapply(aaindex, function(i) {
    vals <- i[["I"]]
    vals <- vals - min(vals, na.rm = TRUE)
    vals/max(vals, na.rm = TRUE)
  }))
  
  colnames(aa_nvals) <- tolower(a(colnames(aa_nvals)))
  aa_nvals
  
  aa_props <- sapply(aaindex, function(i) i[["D"]])
  
  traits <- list(size = c(63, 72, 109, 399),
                 hydroph = c(54, 68, 151, 244),
                 polarity = c(111, 321),
                 alpha = c(3, 41, 253))
  
  # clustering of amino acids 
  
  all_traits_combn <- expand.grid(traits)
  
  all_groups <- apply(all_traits_combn, 1, function(single_trait_combn) {
    #use ward, because Piotr does
    cl <- hclust(dist(t(aa_nvals[single_trait_combn, ])), method = "ward.D2")
    gr <- cutree(cl, k = 4)
    names(gr) <- tolower(names(gr))
    agg_gr <- lapply(unique(gr), function(single_group) names(gr[gr == single_group]))
    names(agg_gr) <- 1L:length(agg_gr)
    agg_gr
  })
  
  
  #interesting encodings 
  
  pos_seqs <- read_uniprot(paste0(pathway, "signal_peptides.txt"), ft_names = "signal")
  
  int_enc <- as.numeric(rownames(p1_dat[p1_dat[, "encoding"] != "", ]))
  group_best <- all_groups[int_enc][[1]]
  #re-arrange group for better comparision
  group_best <- group_best[c(2, 3, 4, 1)]
  names(group_best) <- 1L:4
  
  group_worst <- all_groups[int_enc][[2]]
  
  nhc_borders <- t(sapply(pos_seqs, find_nhc))
  
  deg_pos <- lapply(list(group_best, group_worst), function(single_encoding) {
    single_encoding <- lapply(single_encoding, toupper)
    lapply(1L:length(pos_seqs), function(seq_id) {
      single_seq <- pos_seqs[[seq_id]]
      aa_n <- single_seq[1L:(nhc_borders[seq_id, "start_h"] - 1)]
      aa_h <- single_seq[nhc_borders[seq_id, "start_h"]:(nhc_borders[seq_id, "start_c"] - 1)]
      aa_c <- single_seq[nhc_borders[seq_id, "start_c"]:(nhc_borders[seq_id, "cs"])]
      aa_mature <- single_seq[(nhc_borders[seq_id, "cs"] + 1):length(single_seq)]
      res <- list(aa_n = aa_n,
                  aa_h = aa_h,
                  aa_c = aa_c,
                  aa_m = aa_mature)
      lapply(res, function(i) 
        degenerate(i, single_encoding))
    })
  })
  
  deg_region <- do.call(rbind, lapply(deg_pos, function(single_encoding)
    cbind(region = unlist(lapply(c("n", "h", "c", "m"), function(i) rep(i, 4))), 
          do.call(rbind, lapply(1L:4, function(single_region) {
            res <- factor(unlist(lapply(single_encoding, function(single_seq)
              single_seq[[single_region]])), levels = as.character(1L:4))
            res_tab <- data.frame(table(res))
            cbind(res_tab, prop = res_tab[["Freq"]]/length(res))
          })))))
  
  deg_region <- cbind(enc = unlist(lapply(c("best", "worst"), function(i) rep(i, 16))),
                      deg_region)
  
  colnames(deg_region) <- c("enc", "region", "group", "count", "freq")
  
  deg_region[["region"]] <- factor(deg_region[["region"]], levels = c("n", "h", "c", "m"))
  
  levels(deg_region[["region"]]) <- c("n-region", "h-region", "c-region", "Mature\nprotein")
  levels(deg_region[["group"]]) <- paste0("Group ", levels(deg_region[["group"]]))
  levels(deg_region[["enc"]]) <- c("Best\nsensitivity", "Best\nspecificity")
  
  #save("deg_region", file = "signalHsmm_groups.RData")
  p2 <- ggplot(deg_region, aes(x = region, y = freq, fill = enc, colour = enc)) +
    geom_bar(stat = "identity", position = "dodge") + 
    facet_wrap(~ group, nrow = 1) + 
    scale_y_continuous("Frequency") +
    scale_x_discrete("Region\n") +
    scale_colour_manual("Encoding: ", values = c("red", "blue")) +
    scale_fill_manual("Encoding: ", values = c(adjustcolor("red", 0.25), adjustcolor("blue", 0.25))) + 
    guides(colour = FALSE) +
    my_theme

  #figure comparing encodings ---------------------------
  group_properties <- function(group) {
    res <- do.call(rbind, lapply(1L:length(group), function(subgroup_id)
      melt(data.frame(group = paste0("Group ", as.character(rep(subgroup_id, 4))), critertion = c("size", "hydroph", "polarity", "alpha"),
                      aa_nvals[unlist(all_traits_combn[int_enc[2], ]), group[[subgroup_id]]]))))
    levels(res[["variable"]]) <- toupper(levels(res[["variable"]]))
    res
  }
  
  dat_best <- group_properties(group_best)
  dat_worst <- group_properties(group_worst)
  
  dat_bestworst <- cbind(enc = unlist(lapply(c("best", "worst"), function(i) rep(i, nrow(dat_best)))),
                         rbind(dat_best, dat_worst))
  levels(dat_bestworst[["enc"]]) <- c("Best\nsensitivity", "Best\nspecificity")
  
  
  
  p1 <- ggplot(dat_bestworst, aes(x = critertion, y = value, col = enc, fill = enc)) +
    geom_point(size = 2.5, shape = 21, position = position_dodge(width=0.5)) +
    #geom_text(hjust = -1) +
    facet_grid(~group) +
    scale_x_discrete("Criterion\n", labels = c("size" = "Size","hydroph" = "Hydroph.",
                                               "polarity" = "Polarity","alpha" = expression(paste(alpha, "-helix")))) +
    scale_y_continuous("Value") + 
    scale_colour_manual("Encoding: ", values = c("red", "blue")) +
    scale_fill_manual("Encoding: ", values = c(adjustcolor("red", 0.25), adjustcolor("blue", 0.25))) + 
    my_theme
  
  list(prop_plot = p1, 
       freq_plot = p2, 
       best_sens = group2df(group_best, "The best sensitivity (final) encoding", "tab:best"),
       best_spec = group2df(group_worst, "The best specificity encoding", "tab:worst"))
}

#' Create plot of densities of signal peptide and regions
#'
#' Calculates densities of length of a) whole signal peptides,
#' b) n-regions, c) h-regions d) c-regions.
#'
#' @return a list of two plots. \code{sp_len} is the density of signal peptides
#' lengths. \code{regions} is distribution of regions' lengths.  

plot_reglen <- function() {
  require(reshape2)
  require(signalHsmm)
  require(fitdistrplus)
  require(ggplot2)
  
  fit_dgeom <- function(x) {
    dens <- density(x)
    
    dens_dat <- data.frame(xd = dens[["x"]], yd = dens[["y"]])
    dgeom_dat <- data.frame(x = seq(min(x), max(x))) %>%
      mutate(y = dgeom(x, fitdist(x, "geom")[["estimate"]]))
    
    list(dens_dat = dens_dat,
         dgeom_dat = dgeom_dat)
  }
  
  pos_seqs <- read_uniprot(paste0(pathway, "signal_peptides.txt"), ft_names = "signal")
  
  nhc_borders <- t(sapply(pos_seqs, find_nhc))
  
  # plot of whole signal peptide length
  whole_sp <- data.frame(value = nhc_borders[, "cs"] - 1)
  whole_sp <- cbind(var = rep("Signal peptide", nrow(whole_sp)), whole_sp)
  whole_sp_dat <- fit_dgeom(whole_sp[["value"]])
  whole_sp_dat[["dgeom_dat"]][["var"]] <- "Signal peptide"
  
  p1 <- ggplot(whole_sp_dat[["dgeom_dat"]], aes(x = x, y = y)) +
    geom_bar(stat = "identity", col = "orange", fill = adjustcolor("orange", 0.25)) +
    geom_area(aes(x = xd, y = yd), whole_sp_dat[["dens_dat"]], col = "blue", fill = adjustcolor("blue", 0.25)) + 
    facet_wrap(~ var) +
    scale_y_continuous("Density") + 
    scale_x_continuous("Length\n") + 
    my_theme
  
  # plot of regions' lengths
  
  nhc_len <- sapply(2L:ncol(nhc_borders), function(i)
    nhc_borders[, i] - nhc_borders[, i - 1])
  colnames(nhc_len) <- c("n", "h", "c")
  mlen <- melt(nhc_len)
  levels(mlen[["Var2"]]) <- paste0(levels(mlen[["Var2"]]), "-region")
  
  regions_dat <- lapply(levels(mlen[["Var2"]]), function(single_region) 
    fit_dgeom(mlen[mlen[["Var2"]] == single_region, "value"]))
  
  dgeom_regions <- do.call(rbind, lapply(1L:length(regions_dat), function(single_region_id)
    data.frame(regions_dat[[single_region_id]][["dgeom_dat"]],
               region = levels(mlen[["Var2"]])[single_region_id])))
  
  dens_regions <- do.call(rbind, lapply(1L:length(regions_dat), function(single_region_id)
    data.frame(regions_dat[[single_region_id]][["dens_dat"]],
               region = levels(mlen[["Var2"]])[single_region_id])))
  
  p2 <- ggplot(dgeom_regions, aes(x = x, y = y)) +
    geom_bar(stat = "identity", col = "orange", fill = adjustcolor("orange", 0.25)) +
    geom_area(aes(x = xd, y = yd), dens_regions, col = "blue", fill = adjustcolor("blue", 0.25)) + 
    facet_wrap(~ region, scales = "free_x") +
    scale_x_continuous("Length\n") + 
    scale_y_continuous("Density") + 
    my_theme
  
  list(sp_len = p1, regions = p2)
}

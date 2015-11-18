#' Create plot of densities of signal peptide and regions
#'
#' Calculates densities of length of a) whole signal peptides,
#' b) n-regions, c) h-regions d) c-regions.
#' Fits geometric distribution.
#'
#' @return a list of two plots. \code{sp_len} is the density of signal peptides
#' lengths. \code{regions} is distribution of regions' lengths.  

plot_reglen <- function() {
  get_density <- function(x) {
    dens <- density(x)
    data.frame(xd = dens[["x"]], yd = dens[["y"]])
  }
  
  fit_geom_density <- function(x) {
    data.frame(x = seq(min(x), max(x))) %>%
      mutate(y = dgeom(x, 1/mean(x)))
  }
  
  fit_dgeom <- function(whole_sp, origins = TRUE) {
    #split between plasmodiae and other
    if(origins) {
      dens_dat <- do.call(rbind, lapply(levels(whole_sp[["seq_origin"]]), function(single_origin)
        data.frame(seq_origin = single_origin, get_density(whole_sp[whole_sp[["seq_origin"]] == single_origin, "value"]))))
      dgeom_dat <- do.call(rbind, lapply(levels(whole_sp[["seq_origin"]]), function(single_origin)
        data.frame(seq_origin = single_origin, fit_geom_density(whole_sp[whole_sp[["seq_origin"]] == single_origin, "value"]))))
    } else {
      dens_dat <- get_density(whole_sp[, "value"])
      dgeom_dat <- fit_geom_density(whole_sp[, "value"])
    }
    
    list(dens_dat = dens_dat,
         dgeom_dat = dgeom_dat)
  }
  
  pos_seqs <- read_uniprot(paste0(pathway, "signal_peptides.txt"), ft_names = "signal")
  plas_seqs <- read_uniprot("./plasmodium_benchmark_data/plas.txt", ft_names = "signal")
  seq_origin <- c(rep("Other", length(pos_seqs)), rep("Plasmodium", length(plas_seqs)))
  
  nhc_borders <- t(sapply(c(pos_seqs, plas_seqs), find_nhc))
  
  # plot of whole signal peptide length
  whole_sp <- data.frame(value = nhc_borders[, "cs"] - 1)
  whole_sp <- cbind(seq_origin = seq_origin, var = rep("Signal peptide", nrow(whole_sp)), whole_sp)
  
  whole_sp_dat <- fit_dgeom(whole_sp)
  
  whole_sp_dat[["dgeom_dat"]][["var"]] <- "Signal peptide"
  
  p1 <- ggplot(whole_sp_dat[["dgeom_dat"]], aes(x = x, y = y)) +
    geom_bar(stat = "identity", col = "orange", fill = adjustcolor("orange", 0.25)) +
    geom_area(aes(x = xd, y = yd), whole_sp_dat[["dens_dat"]], col = "blue", fill = adjustcolor("blue", 0.25)) + 
    facet_grid(seq_origin ~ var) +
    scale_y_continuous("Density") + 
    scale_x_continuous("Length\n") + 
    my_theme
  
  # plot of regions' lengths
  
  nhc_len <- sapply(2L:ncol(nhc_borders), function(i)
    nhc_borders[, i] - nhc_borders[, i - 1])
  colnames(nhc_len) <- c("n", "h", "c")
  
  nhc_len <- data.frame(nhc_len, seq_origin = seq_origin)
  mlen <- melt(nhc_len)
  colnames(mlen)[2] <- "Var2"
  levels(mlen[["Var2"]]) <- paste0(levels(mlen[["Var2"]]), "-region")
  
  regions_dat <- lapply(levels(mlen[["Var2"]]), function(single_region) 
    fit_dgeom(mlen[mlen[["Var2"]] == single_region, ]))
  
  dgeom_regions <- do.call(rbind, lapply(1L:length(regions_dat), function(single_region_id)
    data.frame(regions_dat[[single_region_id]][["dgeom_dat"]],
               region = levels(mlen[["Var2"]])[single_region_id])))
  
  dens_regions <- do.call(rbind, lapply(1L:length(regions_dat), function(single_region_id)
    data.frame(regions_dat[[single_region_id]][["dens_dat"]],
               region = levels(mlen[["Var2"]])[single_region_id])))
  
  p2 <- ggplot(dgeom_regions, aes(x = x, y = y)) +
    geom_bar(stat = "identity", col = "orange", fill = adjustcolor("orange", 0.25)) +
    geom_area(aes(x = xd, y = yd), dens_regions, col = "blue", fill = adjustcolor("blue", 0.25)) + 
    facet_grid(seq_origin ~ region, scales = "free_x") +
    scale_x_continuous("Length\n") + 
    scale_y_continuous("Density") + 
    my_theme
  
  list(sp_len = p1, regions = p2)
}

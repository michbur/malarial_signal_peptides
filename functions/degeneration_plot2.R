freq_nondeg <- read.csv2("./frequency_analysis/freq_nondeg.csv")[, -1]
freq_deg <- read.csv2("./frequency_analysis/freq_deg.csv")[, -1]

do_pca <- function(x) 
  x %>% 
  select(-type, -taxon) %>% 
  prcomp(center = TRUE, scale = TRUE) %>% 
  getElement("x") %>% 
  data.frame() %>% 
  select(1, 2) %>% 
  cbind(select(x, type, taxon), .) %>% 
  mutate(type_nice = factor(type, labels = c("mature protein", "signal peptide")),
         taxon_nice = factor(taxon, labels = c("other", "Plasmodium"))) %>% 
  mutate(both = paste0(type_nice, " (", taxon_nice, ")"))

dat_deg <- do_pca(freq_deg) 
dat_nondeg <- do_pca(freq_nondeg)

# ggplot(dat_deg, aes(x = PC1, y = PC2, color = both, shape = both)) + 
#   geom_point(size = 6) 
# 
# ggplot(dat_nondeg, aes(x = PC1, y = PC2, color = both, shape = both)) + 
#   geom_point(size = 6) 

plot_pca <- function(x)
  ggplot(x, aes(x = PC1, y = PC2, color = both, shape = both, fill = both, linetype = both)) + 
  geom_density_2d(color = "black", contour = TRUE) +
  #geom_point() +
  stat_density2d(aes(fill=both,alpha=..level..), color = "black", 
                 contour = TRUE, geom="polygon") +
  scale_linetype_discrete("") +
  scale_fill_discrete("") +
  scale_alpha_continuous(range = c(0.35, 1)) +
  guides(alpha = FALSE) +
  my_theme

p_nondeg <- plot_pca(dat_nondeg) + guides(fill = FALSE, linetype = FALSE) 
p_deg <- plot_pca(dat_deg)

grid.arrange(textGrob("A", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), p_nondeg, 
             textGrob("B", x = 0.75, y = 0.9, gp=gpar(fontsize=22)), p_deg, 
             nrow = 2, ncol = 2, widths = c(0.05, 0.95))

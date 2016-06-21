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
  mutate(both = paste0(type, taxon))

dat_deg <- do_pca(freq_deg) 
dat_nondeg <- do_pca(freq_nondeg)

# ggplot(dat_deg, aes(x = PC1, y = PC2, color = both, shape = both)) + 
#   geom_point(size = 6) 
# 
# ggplot(dat_nondeg, aes(x = PC1, y = PC2, color = both, shape = both)) + 
#   geom_point(size = 6) 

ggplot(dat_deg, aes(x = PC1, y = PC2, color = both, fill = both)) + 
  geom_density_2d() +
  ggtitle("Reduced amino acid alphabet") +
  my_theme

stat_density_2d(aes(x = PC1, y = PC2, color = both, shape = both), dat_nondeg) %>% 
  str

ggplot(dat_nondeg, aes(x = PC1, y = PC2, color = both, shape = both)) + 
  geom_density_2d(contour = TRUE) +
  ggtitle("Full amino acid alphabet") +
  my_theme

setwd("C:/Users/weiga/Dropbox/Bbsl_2023/")
#setwd("Dropbox/Borreliella-U19-Assemblies/clonal-frame-analysis/")
library(tidyverse)

x <- read_tsv("data/Bbss_ldhat_mean_rec.txt")
x %>% ggplot(aes(x=contig, y = mean_rec)) +
  geom_boxplot() + 
  geom_jitter(shape =1, alpha = 0.1) + 
  scale_y_log10() +
  theme_bw()

x %>% group_by(contig) %>% 
  summarise(mean(mean_rec), sd(mean_rec), median(mean_rec))

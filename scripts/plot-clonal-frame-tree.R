setwd("C:/Users/weiga/Dropbox/Borreliella-U19-Assemblies/clonal-frame-analysis/")
#setwd("Dropbox/Borreliella-U19-Assemblies/clonal-frame-analysis/")
library(ape)
#library(genoPlotR)
library(tidyverse)
library(ggtree)
library(readxl)
library(ggrepel)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ggtree")

#############
# Bbss, reconstituted WGS alignments

dd <- read_xlsx("../RicksList-BioSample-v2.xlsx", sheet = "genomes" , col_types = "text") 

plot.ss <- function(prefix) {
  treefile = paste("wgs_", prefix, ".labelled_tree.newick", sep = "")
  recfile = paste("wgs_", prefix, ".importation_status.txt", sep = "")
  
  tree <- treeio::read.tree(treefile)
  
  # recom tracks: remove inodes; map to genome pos
  rec.wide <- read_tsv(recfile)
  
  tips <- tibble(id=tree$tip.label)
  tips <- tips %>% left_join(dd, c("id" = "isolate")) # direct reading didn't work
  col <- c( "green", "red")
  
  p <- ggtree(tree)  %<+% tips +
    geom_tippoint(aes(color = geo)) +
    geom_tiplab(align = TRUE, aes(color = geo), size = 3) +
    scale_color_manual(values = col) 
#    xlim(0,0.01)
  
  facet_plot(p + xlim_tree(0.004), panel="recombination", data = rec.wide, 
             geom = geom_segment, mapping = aes(x = Beg, xend = End, yend = y), 
             color='blue', size = 3) +
    theme_tree2(legend.position=c(.05, .85))
}

plot.ss("cp26") # need to change the function: wgs2_
plot.ss("chr") # wgs_

########################################
# Histogram of recombination tract lengths
glen.chr <- 910724
glen.cp26 <-26498
glen.lp54 <- 53657
glen.lp17 <- 16821
plot.track <- function(prefix, genome_length) {
  recfile <- paste("wgs2_", prefix, ".importation_status.txt", sep = "")
  itv2 <- read_tsv(recfile)
  tlen <- itv2$End-itv2$Beg
  # Identify ones that straddle the original, reapply wrap-around
  wh <- which(itv2$End==genome_length)
  for(i in wh) {
    if(any(itv2$Beg[itv2$Node==itv2$Node[i]]==1)) {
      wh2 = which(itv2$Beg[itv2$Node==itv2$Node[i]]==1)
      tlen[i] = tlen[i]+tlen[itv2$Node==itv2$Node[i]][wh2]
      tlen[itv2$Node==itv2$Node[i]][wh2] = NA
    }
  }

  par(mfrow = c(1,3))
  hist(tlen,100,col="orange3",prob=T)
  hist(log10(tlen),100,col="orange3",prob=T)
  plot.ecdf(log10(tlen),col="orange3")
  par(mfrow=c(1,1))
}

# plot density
itv.chr <- read_tsv("wgs_chr.importation_status.txt")
itv.cp26 <- read_tsv("wgs2_cp26.importation_status.txt")
itv.chr <- itv.chr %>%  mutate(rep = "chromosome")
itv.cp26 <- itv.cp26 %>%  mutate(rep = "cp26")
itv <- bind_rows(itv.chr, itv.cp26)
itv <- itv %>% mutate(track.len = End - Beg + 1)

itv.mean <- itv %>% 
  #  filter(!str_detect(Node, "NODE")) %>% 
  group_by(rep) %>% 
  #  count()
  summarise(med = median(track.len))

itv %>% ggplot(aes(x = track.len)) +
  geom_histogram(binwidth = 50) +
  geom_vline(data = itv.mean, aes(xintercept = med), color = "red", linetype = 2) +
#  scale_y_log10() +
  facet_wrap(~rep, nrow = 2) +
  theme_bw()


########
# rec parameters
dd <- read_xlsx("../RicksList-BioSample-v2.xlsx", sheet = "genomes" , col_types = "text") 

# check alignment quality
x1 <- read_tsv("ss-lp17.win", col_names = c("start", "end", "diff"))
x1 <- x1 %>% mutate(rep = "lp17")

x2 <- read_tsv("ss-cp26.win", col_names = c("start", "end", "diff"))
x2 <- x2 %>% mutate(rep = "cp26")

x3 <- read_tsv("ss-lp54.win", col_names = c("start", "end", "diff"))
x3 <- x3 %>% mutate(rep = "lp54")

x <- bind_rows(x1, x2, x3)
x %>% ggplot(aes(x=start, y=diff)) +
  geom_line()  +
  facet_wrap(~rep, ncol = 1, scales = "free")

# full set of genomes
x1 <- read_tsv("full-lp17.win", col_names = c("start", "end", "diff"))
x1 <- x1 %>% mutate(rep = "lp17")

x2 <- read_tsv("full-cp26.win", col_names = c("start", "end", "diff"))
x2 <- x2 %>% mutate(rep = "cp26")

x3 <- read_tsv("full-lp54.win", col_names = c("start", "end", "diff"))
x3 <- x3 %>% mutate(rep = "lp54")

x <- bind_rows(x1, x2, x3)
x %>% ggplot(aes(x=start, y=diff)) +
  geom_line()  +
  facet_wrap(~rep, ncol = 1, scales = "free")


################
# rec parameters
par <- read_xlsx("../RicksList-BioSample-v2.xlsx", sheet = "ClonalFrame", na = 'NA')
par <- par %>% mutate(sd = sqrt(sd))

par %>% filter(par == 'R/theta' & run == "ss") %>% 
  ggplot(aes(x = replicon, y = mean)) + 
  geom_point(shape = 1, aes(color = type)) +
  theme_bw()

par %>% filter(par == 'r/m' & run == "ss") %>% 
  ggplot(aes(x = replicon, y = mean)) + 
  geom_point(shape = 1, aes(color = type)) +
  theme_bw()


par.delta <- par %>% filter(par == '1/delta' & run == "ss") %>% 
  mutate(track.len = 1/mean) 

par.delta %>% ggplot(aes(x = replicon, y = track.len)) + 
  geom_point(shape = 1, aes(color = type)) +
  theme_bw()
  

par.rt <- par %>% filter(par == 'R/theta') %>% 
  pivot_wider(names_from = "type", values_from = c("mean", "sd"))

par.delta <- par %>% filter(par == '1/delta') %>% 
  pivot_wider(names_from = "run", values_from = c("mean", "sd")) %>% 
  mutate(l.full = 1/mean_full, l.ss = 1/mean_ss)

par.wide %>% ggplot(aes(x = mean_full, y = mean_ss)) +
  geom_pointrange(aes(ymin = mean_ss - sd_ss, ymax = mean_ss + sd_ss)) +
  geom_pointrange(aes(xmin = mean_full - sd_full, xmax = mean_full + sd_full)) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  geom_hline(yintercept = 1) +
  geom_vline(xintercept = 1) +
  geom_point() +
  geom_text_repel(aes(label = replicon)) +
#  scale_x_log10() +
#  scale_y_log10() +
  theme_bw()

par %>% filter(par == 'R/theta') %>% 
  ggplot(aes(x = replicon, y = mean, color = run)) +
  geom_errorbar(aes(ymin = mean -sd, ymax = mean + sd)) +
  geom_hline(yintercept = 1, linetype = 2) +
  theme_bw() +
  geom_point() +
  coord_flip() +
  ylab("Rho/Theta") +
  xlab(label = "") +
  facet_wrap(~type, nrow = 2) +
  theme(legend.position = "bottom")

#############
# iqtrees
iq.chr <- treeio::read.newick("chr-rooted-boot80.dnd")
iq.cp26 <- treeio::read.newick("cp26-rooted-boot80.dnd")
iq.lp54 <- treeio::read.newick("lp17-mid-boot80.dnd")
iq.lp17 <- treeio::read.newick("lp54-rooted-boot80.dnd")

tips <- tibble(id=iq.chr$tip.label)
tips <- tips %>% left_join(dd, c("id" = "isolate")) # direct reading didn't work
col <- c( "green", "blue", "red", "orange")

p.chr <- ggtree(iq.chr)  %<+% tips +
  geom_tippoint(aes(color = geo)) +
  geom_tiplab(align = TRUE, aes(color = geo), size = 3) +
  scale_color_manual(values = col) +
  theme_tree2(legend.position = c(0.05, 0.85)) 

tips <- tibble(id=iq.cp26$tip.label)
tips <- tips %>% left_join(dd, c("id" = "isolate")) # direct reading didn't work
col <- c( "green", "blue", "red", "orange")
p.cp26 <- ggtree(iq.cp26)  %<+% tips +
  geom_tippoint(aes(color = geo)) +
  geom_tiplab(align = TRUE, aes(color = geo), size = 3) +
  scale_color_manual(values = col) +
  theme_tree2(legend.position = c(0.05, 0.85)) 

tips <- tibble(id=iq.lp54$tip.label)
tips <- tips %>% left_join(dd, c("id" = "isolate")) # direct reading didn't work
col <- c( "green", "blue", "red", "orange")
p.lp54 <- ggtree(iq.lp54)  %<+% tips +
  geom_tippoint(aes(color = geo)) +
  geom_tiplab(align = TRUE, aes(color = geo), size = 3) +
  scale_color_manual(values = col) +
  theme_tree2(legend.position = c(0.05, 0.85)) 

tips <- tibble(id=iq.lp17$tip.label)
tips <- tips %>% left_join(dd, c("id" = "isolate")) # direct reading didn't work
col <- c( "green", "blue", "red", "orange")
p.lp17 <- ggtree(iq.lp17)  %<+% tips +
  geom_tippoint(aes(color = geo)) +
  geom_tiplab(align = TRUE, aes(color = geo), size = 3) +
  scale_color_manual(values = col) +
  theme_tree2(legend.position = c(0.05, 0.85)) 

plot_list(p.chr, p.cp26, p.lp54, p.lp17, ncol=2, tag_levels="A")   

# traditional
par(mfrow = c(2,2))
plot(iq.chr, type = 'fan', show.tip.label = T, no.margin = F, main = "chromosome")
add.scale.bar(x=0, y =70)
plot(iq.cp26, type = 'fan', show.tip.label = T, no.margin = F, main = "cp26")
add.scale.bar(x=0, y = 70)
plot(iq.lp54, type = 'fan', show.tip.label = T, no.margin = F, main = "lp54")
add.scale.bar(x =0, y =70)
plot(iq.lp17, type = 'fan', show.tip.label = T, no.margin = F, main = "lp17")
add.scale.bar(x=0, y=70)
par(mfrow=c(1,1))

###############
# plot panels
# chr.length = 910724
plot.rec <- function(prefix) {
  treefile = paste("ss-", prefix, ".labelled_tree.newick", sep = "")
  recfile = paste("ss-", prefix, ".importation_status.txt", sep = "")
  posfile = paste("snp5-", prefix, ".pos", sep = "")
  
  tree <- treeio::read.tree(treefile)
#plot(tree, show.tip.label = F)
#ggtree(tree)  +  
#  geom_treescale() +
#  geom_tiplab(align = T) 
#+
#  xlim(0,0.01)

# recom tracks: remove inodes; map to genome pos
  itv2 = read_tsv(recfile)
  rec.long <- itv2 %>% 
    mutate(seg = 1:nrow(itv2) ) %>% 
    filter(!str_detect(Node, "^NODE")) %>% 
    pivot_longer(2:3, names_to = "coord", values_to = "s.pos")
  pos <- read_tsv(posfile, col_names = F)
  pos.df <- tibble(s.pos = 1:nrow(pos), g.pos = pos$X2)
  rec.long <- rec.long %>% left_join(pos.df, "s.pos")
  rec.wide <- rec.long %>% dplyr::select(1,2,3,5) %>% 
  pivot_wider(names_from = "coord", values_from = "g.pos")

  tips <- tibble(id=tree$tip.label)
  tips <- tips %>% left_join(dd, c("id" = "isolate")) # direct reading didn't work
  col <- c( "green", "blue", "red", "orange")
  
  p <- ggtree(tree)  %<+% tips +
    geom_tippoint(aes(color = geo)) +
    geom_tiplab(align = TRUE, aes(color = geo), size = 3) +
    scale_color_manual(values = col)
  
  facet_plot(p, panel="recombination", data = rec.wide, 
             geom = geom_segment, mapping = aes(x = Beg, xend = End, yend = y), 
             color='blue', size = 3) +
    theme_tree2(legend.position=c(.05, .85))
}

########################################
# Histogram of recombination tract lengths
tlen = rec.wide$End-rec.wide$Beg
# Identify ones that straddle the original, reapply wrap-around
if(FALSE) { # no need
  wh = which(itv2$End==genome_length)
  for(i in wh) {
    if(any(itv2$Beg[itv2$Node==itv2$Node[i]]==1)) {
      wh2 = which(itv2$Beg[itv2$Node==itv2$Node[i]]==1)
      tlen[i] = tlen[i]+tlen[itv2$Node==itv2$Node[i]][wh2]
      tlen[itv2$Node==itv2$Node[i]][wh2] = NA
    }
  }
}

par(mfrow = c(1,3))
hist(tlen,100,col="orange3",prob=T)
hist(log10(tlen),100,col="orange3",prob=T)
plot.ecdf(log10(tlen),col="orange3")
par(mfrow=c(1,1))

#######################
# https://yulab-smu.top/treedata-book/chapter7.html
# panel plot


#tree$tip.label[56:57] <- c("Okinawa-CW61","Okinawa-CW62")
dd <- read_xlsx("../RicksList-BioSample-v2.xlsx", sheet = "genomes" , col_types = "text") 
tips <- tibble(id=tree$tip.label)
tips <- tips %>% left_join(dd, c("id" = "isolate")) # direct reading didn't work
col <- c( "green", "blue", "red", "orange")

p <- ggtree(tree)  %<+% tips +
  geom_tippoint(aes(color = geo)) +
  geom_tiplab(align = TRUE, aes(color = geo), size = 3) +
  scale_color_manual(values = col)

facet_plot(p, panel="recombination", data = rec.wide, 
           geom = geom_segment, mapping = aes(x = Beg, xend = End, yend = y), 
           color='blue', size = 3) +
  theme_tree2(legend.position=c(.05, .85))

# BB_0081 through BB_0084 is at 76451 - 80773
facet_plot(p, panel="recombination", data = rec.wide %>% filter(Beg > 75000 & End < 1e5), 
           geom = geom_segment, mapping = aes(x = Beg, xend = End, yend = y), 
           color='blue', size = 3) +
  theme_tree2(legend.position=c(.05, .85))


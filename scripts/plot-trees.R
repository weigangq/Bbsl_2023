library(tidyverse)

#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("ggtree")

library(ggtree)
#library(ape)
library(readxl)
#library(ggraph)
#library(phytools)

#setwd("Dropbox/Bbsl_2023/")
setwd("c:/Users/weiga/Dropbox/Bbsl_2023/")
# from Lily:
## The codes belowe are used to plot a tree with 2 colored traits and separate legends

###########
# PF54 tree
#############
tree <- treeio::read.tree("pf54-data/pf54-nr.dnd4")
dd <- read_xlsx("../Borreliella-U19-Assemblies/RicksList-BioSample-v2.xlsx", sheet = "pf54-nr", col_types = "text")
tips <- tibble(id=tree$tip.label)
tips <- tips %>% left_join(dd, c("id" = "node")) # direct reading didn't work
col <- c("blue", "cyan", "green", "red", "pink","orange")

p1 <- ggtree(tree) %<+% tips +
  geom_treescale(y=-0.9) + 
  geom_tippoint(aes(color=Geo), size=3) +
  geom_tiplab(align = TRUE) + 
  scale_color_manual(values = col) +
  theme(legend.position = c(0.05, 0.85)) +
  xlim(0,4)

tips2 <- tips %>% select(1,3:25)

tips.long <- tips2 %>% pivot_longer(2:24, names_to = "species", values_to = "strain", values_drop_na = T)

sp.order <- tibble(species = colnames(tips)[3:25], pos = seq(1:23))

tips.long <- tips.long %>% 
  left_join(sp.order, "species")

p1 + geom_facet(panel = "Strains", data = tips.long, geom = geom_text, aes(x = (pos-1)/6, label=strain), hjust = "left") +
  geom_facet(panel = "Strains", data = tips.long, geom = geom_point, aes(x = (pos-1)/6))

############
# species tree
##########

t <- read.tree("chr-snp-23taxa.dnd3")
plot(t, direction = "down", align.tip.label = T, y.lim = 0.5, no.margin = T)

#t <- treeio::read.nhx("species-tree.nhx")
t <- treeio::read.tree("chr_snps-23species.dnd3")
ggtree(t) + 
#  layout_dendrogram() +
#  geom_treescale() +
  geom_tiplab(align = T) +
 xlim(0,1)

# upside down species tree
ggtree(t, branch.length = "none") + 
  layout_dendrogram() +
  geom_tiplab() 
###########
# DbpA tree
#############
tree <- treeio::read.tree("DbpA-data/dbpA-SRC.dnd2")
p1 <- ggtree(tree) +
  geom_treescale(y=-0.9) + 
#  geom_tippoint(aes(color=Geo), size=3) +
  #  scale_color_manual(values = c("#B53CAD", "#A0A0A0")) +
#  geom_tiplab(align = TRUE) + 
#  scale_color_manual(values = col) +
#  theme(legend.position = c(0.05, 0.85)) +
#  xlim(0,2)

dd <- read_xlsx("../Borreliella-U19-Assemblies/RicksList-BioSample-v2.xlsx", sheet = "", col_types = "text")
tips <- tibble(id=tree$tip.label)
tips <- tips %>% left_join(dd, c("id" = "label")) # direct reading didn't work
col <- c("blue", "cyan", "green","pink","red", "orange")

p1 <- ggtree(tree) %<+% tips +
  geom_treescale(y=-0.9) + 
  geom_tippoint(aes(color=Geo), size=3) +
  #  scale_color_manual(values = c("#B53CAD", "#A0A0A0")) +
  geom_tiplab(align = TRUE) + 
  scale_color_manual(values = col) +
  theme(legend.position = c(0.05, 0.85)) +
  xlim(0,2)

tips2 <- tips %>% select(1,4:26)

tips.long <- tips2 %>% pivot_longer(2:24, names_to = "species", values_to = "strain", values_drop_na = T)

sp.order <- tibble(species = colnames(tips)[4:26], pos = seq(1:23))

tips.long <- tips.long %>% 
  left_join(sp.order, "species")

p1 + geom_facet(panel = "Strains", data = tips.long, geom = geom_text, aes(x = (pos-1)/23, label=strain), hjust = "left") +
  geom_facet(panel = "Strains", data = tips.long, geom = geom_point, aes(x = (pos-1)/23))


#### End DbpA tree ##########


###########
# OspC tree
#############

# 1. Read nr tree (cut = 0.05 tree distance)
tree <- treeio::read.tree("OspC-data/oc-SRC.dnd3")
#p1 <- ggtree(tree, layout = "circular") +
#  geom_treescale() + 
#  geom_tippoint(aes(color=Geo), size=3) +
#  geom_tiplab(align = TRUE) 
  

# tree <- treeio::read.tree("OspC-data/OspC-nr.dnd3")
# 2. Read genome meta-data
dd <- read_xlsx("../Borreliella-U19-Assemblies/RicksList-BioSample-v2.xlsx", sheet = "genomes", col_types = "text")

# 3. Read nr sets for OTUs
trim <- read_tsv("OspC-data/oc-SRC.matrix", col_types = paste0(rep("c",25), collapse = ''))

# read species tree:
sp_tree <- treeio::read.newick("chr_snps-23species.dnd3")
# 4. join tables
tips <- tibble(id=tree$tip.label)
tips <- tips %>% 
  mutate(isolate = str_remove(id, "^OspC_")) %>% 
  mutate(strain = str_remove(isolate, "#[123]$") )
tips <- tips %>% left_join(dd, c("strain" = "isolate"))
# make wide table

col <- c("green", "black", "blue", "red","orange")

# I spent the whole Saturday figuring this out: CAN't USE ALL COLUMN. REMOVE!!!!
# It ruined my weekend, literally
tips2 <- tips %>% select(id, geo)

p1 <- ggtree(tree) %<+% tips2 +
  geom_treescale(y=-0.9) + 
  geom_tippoint(aes(color=geo), size=3) +
#  scale_color_manual(values = c("#B53CAD", "#A0A0A0")) +
#  geom_tiplab(align = TRUE, aes(color = new_old)) + 
  scale_color_manual(values = col) +
  theme(legend.position = c(0.05, 0.85)) 
#  xlim(0,1)

#p2 <- p1 + geom_treescale(y=-0.9) + geom_tippoint(aes(color=new_old), size=2) +
#  scale_color_manual(values = c("#B53CAD", "#A0A0A0")) 

#tips2 <- tips %>% select(1,10:32)

tips.long <- trim %>% pivot_longer(2:24, names_to = "species", values_to = "strain", values_drop_na = T)

#tips.long <- tips.long %>% 
#  mutate(x = 1)

sp.order <- tibble(species = sp_tree$tip.label, pos = seq(1:23))
tips.long <- tips.long %>% 
  left_join(sp.order, "species")

tips.long <- tips.long %>% rename(id = otu)
#tips.long2 <- tips.long %>% select(1,3,4)
# assuming y are OTU positions
p1 + 
  geom_facet(panel = "Strains", data = tips.long, geom = geom_text, mapping = aes(x = (pos-1)/25, label = strain), hjust = "left") +
  geom_facet(panel = "Strains", data = tips.long, geom = geom_point, mapping = aes(x = (pos-1)/25)) 

#### End OspC tree ##########


#gheatmap(p1, tips2, offset=0.3,,
#         colnames_angle=-45, hjust=-0.1, width = 0.2)  +
#  scale_fill_grey()
#  scale_fill_manual(values = c("#FF9999", "#FFFF99" , "#CCFF99", "#99FFFF", "#9999FF",  "#00CC66", "#FF99CC", "#FF3333"))
#p <- ggtree(tree) %<+% tips +
#  xlim(0, 1) 

#gheatmap(p, dd, width = 0.4) +
#  scale_fill_manual(breaks=c("0", "1"), values=c("steelblue", "firebrick"))

p +  
  geom_tippoint(aes(color = Geo), size = 3)
  


p <- ggtree(tree) + 
  xlim(0,1) + 
  geom_tiplab(align = T, offset = 0.05) +
  geom_treescale()

col <- c("blue","green","yellow","red","orange")

p1 <- p %<+% tips + geom_tippoint(aes(color=Geo), size=3) +
  scale_color_manual(values = col)

p1  
#  geom_tiplab(aes(color=Geo), size=3, hjust=-0.1)

#p2 <- p1 %<+% tips + geom_tiplab(aes(color=Geo), size=3, hjust=-0.1)
#p12 <- p %<+% tips + geom_tippoint(aes(color=Geo), size=3) +
  

#col <- c("#D95F02","#7570B3","#FF00FF","#FF6600","#0066FF","#1B9E77","#5F5F5F")



#p1x <- p1 + scale_color_manual(values=c("#D95F02","#7570B3","#1B9E77"))
#p2x <- p2 + scale_color_manual(values=c("#FF00FF","#FF6600","#0066FF", "#5F5F5F"))
p12x <- p12 + scale_color_manual(values = col)
p12x

leg1 <- get_legend(p1x)
leg2 <- get_legend(p2x)

pp <- p12x + theme(legend.position="none")
ggdraw(plot_grid(pp, plot_grid(leg1, leg2, nrow=2),
                 rel_widths=c(1, 0.2)))

################################
genome_tree <- read.newick("OspC-data/bb.dnd")
#genome_tree$edge.length <- genome_tree$edge.length * 34
plot(genome_tree)
oc_tree <- read.newick("OspC-data/oc-ss-multi.dnd")
plot(oc_tree)

#t1 <- as_tibble(genome_tree$tip.label)
#t2 <- as_tibble(oc_tree$tip.label)
#asso <- cbind(t1, t2)
otus <- intersect(genome_tree$tip.label, oc_tree$tip.label)

asso <- cbind(otus, otus)

# Ref: http://blog.phytools.org/2015/07/new-method-for-co-phylogenetic-plotting.html
source("http://www.phytools.org/cophylo/v0.1/cophylo.R")
obj <- cophylo(genome_tree, oc_tree, rotate = F, assoc = asso)
plot(obj)


#cophyloplot(genome_tree, oc_tree, assoc = asso, length.line = 4, space = 28, gap = 3, type = "cladogram")

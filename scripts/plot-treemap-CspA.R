library(tidyverse)
library(ggtree)
library(readxl)

# n= 75 genomes have complete PFam54 array
#setwd("Dropbox/Bbsl_2023/")
setwd("c:/Users/weiga/Dropbox/Bbsl_2023/")

# 1. read trimmed gene tree
tree <- treeio::read.tree("pf54-data/cspA-trim-0.3.dnd2")

# 2. Read genome meta-data
dd <- read_xlsx("../Borreliella-U19-Assemblies/RicksList-BioSample-v2.xlsx", sheet = "genomes", col_types = "text")

# 3. Read nr sets for OTUs
trim <- read_tsv("pf54-data/cspA.treemap", col_types = paste0(rep("c",25), collapse = ''))

# 4. read species tree:
sp_tree <- treeio::read.newick("data/chr_snps-23species.dnd3")

# 5. join tables
tips <- tibble(id=tree$tip.label)
tips <- tips %>% 
  mutate(strain = str_remove(id, "^.+\\|") )
tips <- tips %>% left_join(dd, c("strain" = "isolate"))
# make wide table

col <- c("green", "blue", "red", "orange")

# I spent the whole Saturday figuring this out: CAN't USE ALL COLUMN. REMOVE!!!!
# It ruined my weekend, literally
tips2 <- tips %>% select(id, geo)

#p1 <- ggtree(tree) %<+% tips2 +
#  geom_treescale(y=-0.9) + 
#  geom_tippoint(aes(color=geo), size=3) +
  #  scale_color_manual(values = c("#B53CAD", "#A0A0A0")) +
  #  geom_tiplab(align = TRUE, aes(color = new_old)) + 
#  scale_color_manual(values = col) +
#  theme(legend.position = c(0.05, 0.85)) 
#  xlim(0,1)

tips.long <- trim %>% pivot_longer(2:24, names_to = "species", values_to = "strain", values_drop_na = T)

sp.order <- tibble(species = sp_tree$tip.label, pos = seq(1:23))
tips.long <- tips.long %>% 
  left_join(sp.order, "species")

tips.long <- tips.long %>% rename(id = otu)

library(stringr)
tips.long <- tips.long %>% mutate(strain_ct = str_count(strain, '[\\w-]+'))
str_ct <- tips.long %>% group_by(otu) %>% 
  summarise(num_otus = sum(strain_ct)) %>% 
  rename(id = otu)

tips3 <- tips %>% select(id)  %>% 
  left_join(str_ct, "id") %>% 
  mutate(tip.txt = paste(id, "__", num_otus, sep = ''))

tree$tip.label <- tips3$tip.txt
tips.long2 <- tips.long %>% left_join(tips3, c("otu" = "id"))
tips.long2 <- tips.long2 %>% relocate(tip.txt, .before = otu)

tips4 <- tips2 %>% left_join(tips3, "id") %>% 
  relocate(tip.txt, .before = id)

p1 <- ggtree(tree) %<+% tips4 +
  geom_treescale(y=-0.9) + 
  geom_tippoint(aes(color=geo), size=3) +
  scale_color_manual(values = col) +
  theme(legend.position = c(0.05, 0.85)) +
  geom_tiplab(size = 2)

p1 +
  geom_facet(panel = "Strains", data = tips.long2, geom = geom_text, mapping = aes(x = (pos-1)/25, label = strain), hjust = "left", size = 2) +
  geom_facet(panel = "Strains", data = tips.long2, geom = geom_point, mapping = aes(x = (pos-1)/25)) 


################################################################################
#Set up
################################################################################

library(spectrolab)
library(ape)

################################################################################
#Load Data
################################################################################

clean_all = readRDS("clean_all.rds")
all_vn = readRDS("all_vn.rds")

####################################
#Hierarchical Clustering: spectra
####################################

d = dist(as.matrix(clean_all))
hc = hclust(d) 

####################################
#Plot: Spectra
####################################

#create phylogram
hcp = as.phylo(hc)

#change labels
meta = meta(clean_all)
species = unique(meta$Species)
cols = setNames(c("red", "blue", "yellow"), species)

plot(hcp, type = "phylogram", show.tip.label = F, cex = 0.5,
     label.offset = 1)
tiplabels(pch = 19, col = cols)
legend("topleft", names(cols), fill = cols, bty = "n")

####################################
#Hierarchical Clustering: vector normalized spectra
####################################

vnd = dist(as.matrix(all_vn))
vnhc = hclust(vnd) 

####################################
#Plot: vector normalized spectra
####################################

#create phylogram
vn_hcp = as.phylo(vnhc)

#change labels
vn_meta = meta(all_vn)
species = unique(vn_meta$Species)
cols = setNames(c("red", "blue", "yellow"), species)

plot(vn_hcp, type = "phylogram", show.tip.label = F, cex = 0.5,
     label.offset = 1)
tiplabels(pch = 19, col = cols)
legend("topleft", names(cols), fill = cols, bty = "n")
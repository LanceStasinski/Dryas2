library(spectrolab)
library(ape)
library(phytools)
library(vegan)
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
# Distance matrix
################################################################################

#spectra
spec_all = readRDS("Clean-up/Clean_spectra/spec_iS50.rds")
spec_all = spec_all[!meta(spec_all)$Species_ID == "DX",]
spec_all = normalize(spec_all)
spec.mean = aggregate(spec_all,
                      by = meta(spec_all)$Name, 
                      mean, try_keep_txt(mean))


spec.mean = as.matrix(spec.mean)

#change row names
meta_new = read.csv("is50_meta.csv", stringsAsFactors = F)
rownames(spec.mean) <- meta_new[, "genetic_name"]


#Genetic distance matrix
tree = read.tree(file = "dryas_phylogeny.treefile")
tree.dist = cophenetic(tree)

#Make sure both matrices have same sample IDs
gene.ID = colnames(tree.dist)
spec.ID = rownames(spec.mean)

remove.gene = setdiff(gene.ID, spec.ID)
remove.spec = setdiff(spec.ID, gene.ID)

tree.dist = tree.dist[!rownames(tree.dist) %in% remove.gene, ]
tree.dist = tree.dist[,!colnames(tree.dist) %in% remove.gene]
spec.mean = spec.mean[!rownames(spec.mean) %in% remove.spec,]

#reorder spec matrix to match order of genetic matrix
spec.mean = spec.mean[row.names(tree.dist),,drop = F]

#calculate distance
spec.dist = as.matrix(dist(spec.mean,
                           method = "manhattan", 
                           diag = T, upper = T))


#log transform spectral distances
spec.dist = log10(spec.dist)

################################################################################
#compare matrices
################################################################################
spec.dist2 = as.dist(spec.dist)
tree.dist2 = as.dist(tree.dist)

mtest = vegan::mantel(spec.dist2, tree.dist2, method = "spearman", 
                      permutations = 9999)
mtest
################################################################################
#plot
################################################################################
library(hexbin)
library(grid)

bin = hexbin(log10(spec.dist2), tree.dist2, xbins = 50)
plot(bin, main = "Genetic distance versus spectral distance",
     xlab = "",
     ylab = "", colramp=BTY)
grid.text("log10(Spectral distance)", .45, .07, gp=gpar(fontsize=16))
grid.text("Genetic distance", .015, .5, rot=90, gp=gpar(fontsize=16))

################################################################################
#heirarchical clustering
################################################################################
library(fpc)
library(cluster)

pam = pamk(spec.dist)
pam #2 clusters recommended

fit = hclust(spec.dist2, method = 'ward.D2')
colors = c("red", "blue")
clus2 = cutree(fit, 2)
plot(as.phylo(fit), type = "fan", tip.color = colors[clus2],
     label.offset = .1, cex = 0.7)

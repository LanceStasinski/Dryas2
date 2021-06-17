library(spectrolab)
library(ape)
library(phytools)
library(vegan)
library(parameters)
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
#Average spectra per individual and phylogeny with same name
################################################################################

#spectra
spec_all = readRDS("Clean-up/Clean_spectra/clean_all_6scans.rds")
spec_all = spec_all[!meta(spec_all)$Species_ID == "DX",]
spec.mean = aggregate(spec_all,
                      by = meta(spec_all)$Name, 
                      mean, try_keep_txt(mean))
spec.mean = as.matrix(spec.mean)


spec.se = aggregate(spec_all,
                      by = meta(spec_all)$Name, 
                      parameters::standard_error, try_keep_txt(mean))
spec.se = as.matrix(spec.se)


#change row names
meta_new = read.csv("mean_meta_new.csv", stringsAsFactors = F)
rownames(meta_new) = meta_new$Name
meta_new = meta_new[row.names(spec.mean),,drop = F]
rownames(spec.mean) <- meta_new[, "genetic_name"]
rownames(spec.se) <- meta_new[, "genetic_name"]


#Phylogeny
tree = read.tree(file = "dryas_phylogeny.treefile")

#make sure both phylogeny and spectra represent same individuals
gene.ID = tree$tip.label
spec.ID = rownames(spec.mean)

remove.gene = setdiff(gene.ID, spec.ID)
remove.spec = setdiff(spec.ID, gene.ID)

tree_clean = drop.tip(tree, remove.gene)
spec.mean = spec.mean[!rownames(spec.mean) == "tmi43_DA",]
spec.se = spec.se[!rownames(spec.se) == "tmi43_DA",]


#reorder spec matrix to match order of genetic matrix
spec.mean = spec.mean[tree_clean$tip.label,,drop = F]
spec.se = spec.se[tree_clean$tip.label,,drop = F]

spec = as_spectra(spec.mean)
se = as_spectra(spec.se)

################################################################################
#Phylogenetic signal
################################################################################

p = c()
l = c()

for (i in seq(400, 2400, 1)) {
  signal = phylosig(tree_clean, spec[,i], method = 'K', test = T, nsim = 1000,
                    se = se[,i])
  lambda = signal$lambda
  l = append(l, lambda)
  p_val = signal$P
  p = append(p, p_val)
}

################################################################################
#plot
################################################################################

par(mfrow = c(2,1))

plot(l, 
     xlab = 'Wavelength (nm)',
     ylab = 'lambda',
     xaxt = 'n')
axis(1,
     at = c(0, 500, 1000, 1500, 2000), 
     labels = c(400, 900, 1400, 1900, 2400))

plot(p,
     xlab = 'Wavelength (nm)',
     ylab = 'p-value',
     xaxt = 'n')
axis(1,
     at = c(0, 500, 1000, 1500, 2000), 
     labels = c(400, 900, 1400, 1900, 2400))





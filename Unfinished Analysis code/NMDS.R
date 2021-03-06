################################################################################
#load package
################################################################################

library(vegan)

################################################################################
#Set up
################################################################################

spectra = readRDS("clean_all.rds")
spectra = resample(spectra, seq(400, 2400, by = 10))
names(spectra) = meta(spectra)$Species
spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)

################################################################################
#NMMDS
################################################################################

mds = metaMDS(spectra.m)
names(mds)

vscores = mds$species
sscores = mds$points

plot(mds)

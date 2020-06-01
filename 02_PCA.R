#PCA
#RDS files for cleaned data can be downloaded from:
#FM_Polar_Studies > Dryas_Spectral_Analyses > Clean-up 
#or generated via 01_clean_up.R
################################################################################
#Set up
################################################################################

library(mixOmics)
library(spectrolab)
library(ggplot2)
library(ggfortify)

################################################################################
#Set working directory to folder containing downloaded rds files
################################################################################

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Clean-up")

################################################################################
#PCA all spectra
################################################################################

#load data
spectra = readRDS("Clean_spectra/clean_all.rds")
spectra = resample(spectra, seq(400, 2400, by = 10))
names(spectra) = meta(spectra)$Species
spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)

#PCA
spectra.pca = prcomp(spectra.m, scale =TRUE)

autoplot(spectra.pca, data = spectra.df, colour = "Species", loadings = TRUE, 
         loadings.colour = 'blue', loadings.label = TRUE, 
         loadings.label.size =3, frame = TRUE, frame.type = 'norm',
         x = 1, y = 3)

################################################################################
#PCA all spectra: vector normalized
################################################################################

#load data
spectra = readRDS("Vector_normalized/all_vn.rds")
spectra = resample(spectra, seq(400, 2400, by = 10))
names(spectra) = meta(spectra)$Species
spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)

#PCA
spectra.pca = prcomp(spectra.m, scale =TRUE)

autoplot(spectra.pca, data = spectra.df, colour = "Species", loadings = TRUE, 
         loadings.colour = 'blue', loadings.label = TRUE, 
         loadings.label.size =3, frame = TRUE, frame.type = 'norm',
         x = 1, y = 3)

################################################################################
#3 sites, all species
################################################################################
spectra = readRDS("Clean_spectra/clean_big3.rds")
spectra = resample(spectra, seq(400, 2400, by = 10))
names(spectra) = meta(spectra)$Species
spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)

#PCA
spectra.pca = prcomp(spectra.m, scale =TRUE)

autoplot(spectra.pca, data = spectra.df, colour = "Species", loadings = TRUE, 
         loadings.colour = 'blue', loadings.label = TRUE, 
         loadings.label.size =3, frame = TRUE, frame.type = 'norm',
         x = 3, y = 4)

################################################################################
#3 sites, no hybrid
################################################################################

spectra = readRDS("Clean_spectra/big3.no_hybrids.rds")
spectra = resample(spectra, seq(400, 2400, by = 10))
names(spectra) = meta(spectra)$Species
spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)

#PCA
spectra.pca = prcomp(spectra.m, scale =TRUE)

autoplot(spectra.pca, data = spectra.df, colour = "Species", loadings = TRUE, 
         loadings.colour = 'blue', loadings.label = TRUE, 
         loadings.label.size =3, frame = TRUE, frame.type = 'norm',
         x = 4, y = )

################################################################################
#PCA octopetala by site
################################################################################

#load data
spectra = readRDS("Clean_spectra/clean_all.rds")
octopetala = spectra[meta(spectra)$Species == "octopetala",]
spectra = octopetala
spectra = resample(spectra, seq(400, 2400, by = 10))
names(spectra) = meta(spectra)$Location
spectra.df = as.data.frame(spectra)
spectra.m = as.matrix(spectra)

#PCA
spectra.pca = prcomp(spectra.m, scale =TRUE)

autoplot(spectra.pca, data = spectra.df, colour = "Location", loadings = TRUE, 
         loadings.colour = 'blue', loadings.label = TRUE, 
         loadings.label.size =3, frame = TRUE, frame.type = 'norm',
         x = 1, y = 3)


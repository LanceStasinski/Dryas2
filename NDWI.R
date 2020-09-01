#NDWI

library(spectrolab)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
# obtain NDWI
################################################################################

#PSR+
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
data = meta(spec_all)

ndwi = (spec_all[,860] - spec_all[,1240])/(spec_all[,860] + spec_all[,1240])

data = cbind(data, ndwi)

boxplot(ndwi ~ Location, data = data, ylab = 'NDWI')

a = aov(ndwi ~ Location, data = data)
summary(a)

#iS50
ispec = readRDS('Clean-up/Clean_spectra/spec_iS50.rds')
spec_all = ispec
spec_m = as.matrix(spec_all)
data = meta(spec_all)
ndwi = (spec_all[,859.99312005504] - spec_all[,1199.90400767939])/(spec_all[,859.99312005504] + spec_all[,1199.90400767939])
data = cbind(data, ndwi)

boxplot(ndwi ~ Location, data = data, ylab = 'NDWI')
a = aov(ndwi ~ Location, data = data)
summary(a)
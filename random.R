library(spectrolab)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
#Data setup !!!!NOTE: use ctrl+f to find a replace the field to be classified!!!
################################################################################

#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")

data = meta(spec_all)

boxplot(LMA~Species_ID, data = data)

wdb = data[data$Location == 'Wickersham Dome B',]

boxplot(LMA~Species_ID, data = wdb)

dx = spec_all[meta(spec_all)$Species_ID == 'DX',]

dev.new(width = 6, height = 8, unit = 'in')
par(mfrow = c(1,2))
boxplot(DA~Species_ID, data = spec_all, 
        names = c('D. alaskensis', 'D. ajanensis', 'Hybrid'), 
        xlab = "Species",
        ylab = 'Dryas alaskensis ancestry')
boxplot(DA~Location, data = dx,
        xlab = "Location",
        ylab = 'Dryas alaskensis ancestry')
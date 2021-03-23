#Plot Spec
library("spectrolab")

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

#vector normalized

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
spec_all = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_wet = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
ala = spec_all[meta(spec_all)$Species_ID == "DA",]
oct = spec_all[meta(spec_all)$Species_ID == "DO",]
hyb = spec_all[meta(spec_all)$Species_ID == "DX",]

es_dry = spec_all[meta(spec_all)$Location == "Eagle Summit",]
es_wet = spec_wet[meta(spec_wet)$Location == "Eagle Summit",]

dry = es_dry[meta(es_dry)$Name == "DryoctESA52",]
wet = es_wet[meta(es_wet)$Name == "DryoctESA52",]

es_wet_oct = es_wet[meta(es_wet)$Species_ID == "DO",]

par(mfrow=c(1,2), mar=c(4,4,4,2), oma = c(4,3,2,2))
plot(wet, lwd = 0.5, lty = 1, col = "grey50", main="Wet", 
     cex.lab = 1.2, ylim = c(0, .04), ylab = "Vector Normalized Reflectance", xlab = "Wavelength(nm)")
plot_regions(wet, regions = default_spec_regions(), add = TRUE)

plot(dry, lwd = 0.5, lty = 1, col = "grey50", main="Dry", 
     cex.lab = 1.2, ylim = c(0, .04), ylab = NA, xlab = "Wavelength(nm)")
plot_regions(dry, regions = default_spec_regions(), add = TRUE)


par(mfrow=c(1,3), mar=c(4,4,4,2), oma = c(4,3,2,2))
plot(mean(ala), lwd = 0.5, lty = 1, col = "grey50", main="Alaskensis", 
     cex.lab = 1.2, ylim = c(0, .04), ylab = "Vector Normalized Reflectance", xlab = NA)
plot_quantile(ala, total_prob = 0.8, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(ala, regions = default_spec_regions(), add = TRUE)

plot(mean(oct), lwd = 0.5, lty = 1, col = "grey50", main="Octopetala", 
     cex.lab = 1.2,  xlab = "Wavelength (nm)", ylim = c(0, .04), ylab = NA)
plot_quantile(oct, total_prob = 0.8, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(oct, regions = default_spec_regions(), add = TRUE)

plot(mean(hyb), lwd = 0.5, lty = 1, col = "grey50", main="Hybrid", xlab = NA,
     ylim = c(0, .04), ylab = NA)
plot_quantile(hyb, total_prob = 0.8, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(hyb, regions = default_spec_regions(), add = TRUE)

mtext("Mean spectra and 80% spectral quantile -  Dry Scans", outer=TRUE,  cex=1, line=-0.5)


#Regular Scans

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
ala = spec_all[meta(spec_all)$Species_ID == "DA",]
oct = spec_all[meta(spec_all)$Species_ID == "DO",]
hyb = spec_all[meta(spec_all)$Species_ID == "DX",]

par(mfrow=c(1,3), mar=c(4,4,4,2), oma = c(4,3,2,2))
plot(mean(ala), lwd = 0.5, lty = 1, col = "grey50", main="Alaskensis", 
     cex.lab = 1.2, ylim = c(0, .85), ylab = "Reflectance", xlab = NA)
plot_quantile(ala, total_prob = 0.95, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(ala, regions = default_spec_regions(), add = TRUE)

plot(mean(oct), lwd = 0.5, lty = 1, col = "grey50", main="Octopetala", 
     cex.lab = 1.2,  xlab = "Wavelength (nm)", ylim = c(0, .85), ylab = NA)
plot_quantile(oct, total_prob = 0.95, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(oct, regions = default_spec_regions(), add = TRUE)

plot(mean(hyb), lwd = 0.5, lty = 1, col = "grey50", main="Hybrid", xlab = NA, ylim = c(0, .85),
    ylab = NA)
plot_quantile(hyb, total_prob = 0.95, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(hyb, regions = default_spec_regions(), add = TRUE)

mtext("Mean spectra and 95% spectral quantile", outer=TRUE,  cex=1, line=-0.5)


#hybrid variation vs parental variation
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
ala = spec_all[meta(spec_all)$Species_ID == "DA",]
oct = spec_all[meta(spec_all)$Species_ID == "DO",]
hyb = spec_all[meta(spec_all)$Species_ID == "DX",]

sd.hyb = sd(hyb)/mean(spec_all)
sd.ala = sd(ala)/mean(spec_all)
sd.oct = sd(oct)/mean(spec_all)

plot(sd.ala, lwd = 2, lty = 1, col = "#00B0F6", 
     main = "Mean reflectance per species", cex.lab = 1.5,
     ylim = range(value(sd.hyb),value(sd.ala), value(sd.oct)),
     ylab = "Reflectance", xlab = "Wavelength (nm)")
plot(sd.hyb, lwd = 2, lty = 1, col = "black", add = T)
plot(sd.oct, lwd = 2, lty = 1, col = "#F8766D", add = T)
legend('bottomright', inset = .02, legend = c("DA", "DO", "DX"),
       col = c("#00B0F6", "#F8766D", "black"), lty = 1, lwd = 2, cex =.8,
       bg = "white")

#Just mean reflectance
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
ala = spec_all[meta(spec_all)$Species_ID == "DA",]
oct = spec_all[meta(spec_all)$Species_ID == "DO",]
hyb = spec_all[meta(spec_all)$Species_ID == "DX",]

dev.new(width = 6, height = 6, unit = 'in')
plot(mean(ala), lwd = 2, lty = 1, col = "#00B0F6", 
     main = "Mean reflectance per species", cex.lab = 1.5, ylim = c(0, .85),
     ylab = "Reflectance", xlab = "Wavelength (nm)")
plot(mean(oct), lwd = 2, lty = 1, col = "#F8766D", add = T)
plot(mean(hyb), lwd = 2, lty = 1, col = "black", add = T)
legend('bottomright', inset = .02, legend = c("DAK", "DAJ", "DX"),
       col = c("#00B0F6", "#F8766D", "black"), lty = 1, lwd = 2, cex =.8,
       bg = "white")

#Populations
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
s.m = as_spectra(as.matrix(spec_all))
meta(s.m) = read.csv('metadata_2.csv', stringsAsFactors = F)
spec_all = s.m


da_et = spec_all[meta(spec_all)$GenePop_ID == "DA_et",]
da_wd = spec_all[meta(spec_all)$GenePop_ID == 'DA_wdb']
et = spec_all[meta(spec_all)$GenePop_ID == "DO_et",]
wd = spec_all[meta(spec_all)$GenePop_ID == "DO_wd",]
mdb = spec_all[meta(spec_all)$GenePop_ID == "DO_mdb",]
bg = spec_all[meta(spec_all)$GenePop_ID == "DO_bgc",]

dev.new(width = 6, height = 8, unit = 'in')
plot(mean(da_et), lwd = 2, lty = 1, col = "#a6cee3", 
     cex.lab = 1.5, ylim = c(0, .85), ylab = "Reflectance", 
     xlab = 'Wavelength (nm)', main = "Mean reflecatnce per population")
plot(mean(da_wd), lwd = 2, lty = 1, col = "#b2df8a", add = T)
plot(mean(et), lwd = 2, lty = 1, col = "#fb9a99", add = T)
plot(mean(wd), lwd = 2, lty = 1, col = "#ff7f00", add = T)
plot(mean(mdb), lwd = 2, lty = 1, col = "#e31a1c", add = T)
plot(mean(bg), lwd = 2, lty = 1, col = "#33a02c", add = T)
legend('bottomright',inset = .02, 
       legend=c('DAK ES-TM', 'DAK WDB', 'DAJ BG', 'DAJ ES-TM', 'DAJ MD', 'DAJ WD'),
       col=c("#a6cee3",'#b2df8a', "#33a02c", "#fb9a99", "#e31a1c", "#ff7f00"), 
       lty=1, lwd = 2, cex=0.8, bg = 'white')





#?
par(mfrow = c(1,1))
plot(mean(ala), lwd = 0.5, lty = 1, col = 'blue', 
     cex.lab = 1.2, ylim = c(0, .85), ylab = "Reflectance", xlab = "Wavelength (nm)")
plot(mean(hyb), lwd = 0.5, lty = 1, col = 'black', add = T)
plot(mean(oct), lwd = 0.5, lty = 1, col = 'red', add = T)
plot_regions(ala, regions = default_spec_regions(), add = TRUE)
legend('bottomright',inset = .02, legend=c("D. alaskensis", "Hybrid", "D. octopetala"),
col=c('blue', 'gray50', 'red'), text.font = c(3,1,3), lty=1, cex=0.8)





plot_quantile(ala, total_prob = 0.6, col = rgb(0, 0, 1, 0.25), border = FALSE, add = TRUE)
plot_quantile(hyb, total_prob = 0.6, col = 'gray95', border = FALSE, add = TRUE)
plot_quantile(ala, total_prob = 0.6, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)

#oct wet vs dry
spec_wet = readRDS("Clean-up/Clean_spectra/clean_all_w.rds")
wet = spec_wet[meta(spec_wet)$Species_ID == "DO",]

spec_dry = readRDS("Clean-up/Clean_spectra/clean_all.rds")
dry = spec_dry[meta(spec_dry)$Species_ID == "DO",]

par(mfrow=c(1,2), mar=c(4,4,4,2), oma = c(4,3,2,2))
plot(mean(wet), lwd = 0.5, lty = 1, col = "grey50", main="Wet", 
     cex.lab = 1.2,  xlab = "Wavelength (nm)", ylim = c(0, .85), ylab = "Reflectance")
plot_quantile(wet, total_prob = 0.8, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(wet, regions = default_spec_regions(), add = TRUE)

plot(mean(dry), lwd = 0.5, lty = 1, col = "grey50", main="Dry", 
     cex.lab = 1.2, ylim = c(0, .85), ylab = "Reflectance", xlab = "Wavelength (nm)")
plot_quantile(dry, total_prob = 0.8, col = rgb(1, 0, 0, 0.25), border = FALSE, add = TRUE)
plot_regions(dry, regions = default_spec_regions(), add = TRUE)


mtext("Octopetala Wet versus Dry Scans", outer=TRUE,  cex=1, line=-0.5)

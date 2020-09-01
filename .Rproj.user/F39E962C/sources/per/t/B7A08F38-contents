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

#Populations
da = spec_all[meta(spec_all)$GenePop_ID == "DA",]
et = spec_all[meta(spec_all)$GenePop_ID == "DO_et",]
wdb = spec_all[meta(spec_all)$GenePop_ID == "DO_wdb",]
mdb = spec_all[meta(spec_all)$GenePop_ID == "DO_mdb",]
wda = spec_all[meta(spec_all)$GenePop_ID == "DO_wda",]
bg = spec_all[meta(spec_all)$GenePop_ID == "DO_bgc",]

dev.new(width = 6, height = 6, unit = 'in')
par(mar=c(5,5,2,2))
plot(mean(da), lwd = 2, lty = 1, col = "red", 
     cex.lab = 1.5, ylim = c(0, .85), ylab = "Reflectance", xlab = 'Wavelength (nm)')
plot_regions(da, regions = default_spec_regions(), add = TRUE)
plot(mean(et), lwd = 2, lty = 1, col = "green3", add = T)
plot(mean(wdb), lwd = 2, lty = 1, col = "darkviolet", add = T)
plot(mean(mdb), lwd = 2, lty = 1, col = "cyan2", add = T)
plot(mean(wda), lwd = 2, lty = 1, col = "brown", add = T)
plot(mean(bg), lwd = 2, lty = 1, col = "orange", add = T)
legend('bottomright',inset = .02, 
       legend=c('DA', 'DO_bg', 'DO_et', 'DO_mdb', 'DO_wda', 'DO_wdb'),
       col=c('red', 'orange', 'green3', 'cyan2', 'brown', 'darkviolet'), 
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

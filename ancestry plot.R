setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
d = meta(spec_all)
w = subset(d, Location == "Wickersham Dome B")
e = subset(d, Location == "Eagle Summit")
t = subset(d, Location == "Twelve Mile")
par(mfrow = c(1,4))
boxplot(w$DA~w$Species_ID, main = 'WDB', xlab = 'Species',
        ylab = 'Proportion of ancestry from DA', cex = 1.5)
boxplot(e$DA~e$Species_ID, main = 'ES', xlab = 'Species', ylab = NA)
boxplot(t$DA~t$Species_ID, main = "TM", xlab = 'Species', ylab = NA)
boxplot(d$DA~d$Species_ID, main = "All Sites", xlab = 'Species', ylab = NA)

#NDWI

library(spectrolab)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")


################################################################################
# obtain NDWI
################################################################################
#PSR+
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
do = spec_all[meta(spec_all)$Species_ID == 'DA',]
data = meta(do)

ndwi = (do[,860] - do[,1240])/(do[,860] + do[,1240])

data = cbind(data, ndwi)

boxplot(ndwi ~ Location, data = data, ylab = 'NDWI')

pairwise.t.test(data$ndwi, data$Location)

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



################################################################################
# Water bands
################################################################################
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
do = spec_all[meta(spec_all)$Species_ID == 'DO',]

bg = do[meta(do)$Location == 'Bison Gulch',]
bg = bg[,400]
bg.df = as.data.frame(bg)
bgn = rep('bg', length(bg.df))
bg.df = cbind(bg.df, bgn)
colnames(bg.df) <- c("reflectance", 'location')

es = do[meta(do)$Location == 'Eagle Summit',]
es = es[,400]
es.df = as.data.frame(es)
esn = rep('es', length(es.df))
es.df = cbind(es.df, esn)
colnames(es.df) <- c("reflectance", 'location')

md = do[meta(do)$Location == 'Murphy Dome B',]
md = md[,400]
md.df = as.data.frame(md)
mdn = rep('md', length(md.df))
md.df = cbind(md.df, mdn)
colnames(md.df) <- c("reflectance", 'location')

tm = do[meta(do)$Location == 'Twelve Mile',]
tm = tm[,400]
tm.df = as.data.frame(tm)
tmn = rep('tm', length(tm.df))
tm.df = cbind(tm.df, tmn)
colnames(tm.df) <- c("reflectance", 'location')

wda = do[meta(do)$Location == 'Wickersham Dome A',]
wda = wda[,400]
wda.df = as.data.frame(wda)
wdan = rep('wda', length(wda.df))
wda.df = cbind(wda.df, wdan)
colnames(wda.df) <- c("reflectance", 'location')

wdb = do[meta(do)$Location == 'Wickersham Dome B',]
wdb = wdb[,400]
wdb.df = as.data.frame(wdb)
wdbn = rep('wdb', length(wdb.df))
wdb.df = cbind(wdb.df, wdbn)
colnames(wdb.df) <- c("reflectance", 'location')

df = Reduce(rbind, list(bg.df, es.df, md.df, tm.df, wda.df, wdb.df))

boxplot(reflectance~location, data = df)

pairwise.t.test(df$reflectance, df$location)



spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
do = spec_all[meta(spec_all)$Species_ID == 'DA',]

es = do[meta(do)$Location == 'Eagle Summit',]
es = es[,1940]
es.df = as.data.frame(es)
esn = rep('es', length(es.df))
es.df = cbind(es.df, esn)
colnames(es.df) <- c("reflectance", 'location')

tm = do[meta(do)$Location == 'Twelve Mile',]
tm = tm[,1940]
tm.df = as.data.frame(tm)
tmn = rep('tm', length(tm.df))
tm.df = cbind(tm.df, tmn)
colnames(tm.df) <- c("reflectance", 'location')

wdb = do[meta(do)$Location == 'Wickersham Dome B',]
wdb = wdb[,1940]
wdb.df = as.data.frame(wdb)
wdbn = rep('wdb', length(wdb.df))
wdb.df = cbind(wdb.df, wdbn)
colnames(wdb.df) <- c("reflectance", 'location')

df = Reduce(rbind, list(es.df, tm.df, wdb.df))

boxplot(reflectance~location, data = df)

p =pairwise.t.test(df$reflectance, df$location)
p$p.value

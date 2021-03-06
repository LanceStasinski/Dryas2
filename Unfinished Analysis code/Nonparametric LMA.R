#LMA

################################################################################
#Load Data
################################################################################
library(spectrolab)
library(dplyr)
library(ggpubr)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0/Clean-up/Clean_spectra")
clean_all = readRDS("clean_all.rds")

################################################################################
#Create LMA objects
################################################################################
#All LMA
lma_all = meta(clean_all)$LMA

#LMA by species
lma_oct = meta(clean_all[meta(clean_all)$Species == "octopetala", ])$LMA
lma_ala = meta(clean_all[meta(clean_all)$Species == "alaskensis", ])$LMA
lma_hyb = meta(clean_all[meta(clean_all)$Species == "hybrid", ])$LMA

#LMA by Location
lma_tm = meta(clean_all[meta(clean_all)$Location == "Twelve Mile", ])$LMA
lma_es = meta(clean_all[meta(clean_all)$Location == "Eagle Summit", ])$LMA
lma_mdb = meta(clean_all[meta(clean_all)$Location == "Murphy Dome B", ])$LMA
lma_bg = meta(clean_all[meta(clean_all)$Location == "Bison Gulch", ])$LMA
lma_wda = meta(clean_all[meta(clean_all)$Location == "Wickersham Dome A", ])$LMA
lma_wdb = meta(clean_all[meta(clean_all)$Location == "Wickersham Dome B", ])$LMA

################################################################################
#Investigate normality of LMA values
#Will be used to determine proper statistical test to run
################################################################################

#Density plots
ggdensity(lma_all, 
          main = "Density plot of All LMA",
          xlab = "LMA") #moderately bell shaped


ggdensity(lma_oct, 
          main = "Density plot of Oct LMA",
          xlab = "LMA") #close to bell shaped


ggdensity(lma_ala, 
          main = "Density plot of Ala LMA",
          xlab = "LMA") #less bell shaped

ggdensity(lma_hyb, 
          main = "Density plot of Hyb LMA",
          xlab = "LMA") #nearly bimodal - might be due to half of the samples
                        #being closely related to ala and the other half oct
#QQ plots
ggqqplot(lma_all) #deviates from normality

ggqqplot(lma_ala) #deviates form normality

ggqqplot(lma_oct) #deviates from normality

ggqqplot(lma_hyb) #deviates form normality

#Shapiro-wilks test
shapiro.test(lma_all) #p = 2.2e-16 --- not normal
shapiro.test(lma_oct) #p = 1.692e-15 --- not normal
shapiro.test(lma_ala) #p = 9.684e-16 --- not normal
shapiro.test(lma_hyb) #p = .008101 --- not normal


################################################################################
#Nonparametric tests
#The test above show that LMA values are not normally distributed. As such, 
#parametric tests will not be reliable evaluations of the data. 
################################################################################

#Wilcox test
all = clean_all[!meta(clean_all)$LMA == "NA",]
octala = all[!meta(all)$Species == "hybrid",]
alahyb = all[!meta(all)$Species == "octopetala",]
octhyb = all[!meta(all)$Species == "alaskensis",]

wilcox.test(LMA ~ Species, data = meta(octala), paired = F)#p = 1.5e-13 =dif pops
wilcox.test(LMA ~ Species, data = meta(alahyb), paired =F)#p = .015 = dif pops
wilcox.test(LMA ~ Species, data = meta(octhyb), paired = F)#p = .5977 = NOT DIF

#According to the wilcox test, the LMA for alaskensis is sig. dif. than the LMA 
#of octopetala and hybrid. However, the LMA for hybrids and octopetala are not 
#significantly different. 


wilcox = function(x,locations){
  lma1 = x[meta(x)$Location == locations,]
  wilcox.test(LMA ~ Location, data = meta(lma1), paired = F)
}




oct = clean_all[meta(clean_all)$Species == "octopetala", ]

es_tm = c("Eagle Summit", "Twelve Mile")
es_bg = c("Eagle Summit", "Bison Gulch")
es_mdb = c("Eagle Summit", "Murphy Dome B")
es_wda = c("Eagle Summit", "Wickersham Dome A")
es_wdb = c("Eagle Summit", "Wickersham Dome B")
tm_bg = c("Twelve Mile", "Bison Gulch")
tm_mdb = c("Twelve Mile", "Murphy Dome B")
tm_wda = c("Twelve Mile", "Wickersham Dome A")
tm_wdb = c("Twelve Mile", "Wickersham Dome B")
bg_mdb = c("Bison Gulch", "Murphy Dome B")
bg_wda = c("Bison Gulch", "Wickersham Dome A")
bg_wdb = c("Bison Gulch", "Wickersham Dome B")
mdb_wda = c("Murphy Dome B", "Wickersham Dome A")
mdb_wdb = c("Murphy Dome B", "Wickersham Dome B")
wda_wdb = c("Wickersham Dome A", "Wickersham Dome B")


o_es_tm = wilcox(oct, locations = es_tm) #p = .1006
o_es_bg = wilcox(oct, locations = es_bg) #p = .0003692
o_es_mdb = wilcox(oct, locations = es_mdb) #p = .001633
o_es_wda = wilcox(oct, locations = es_wda) #p = 1.294e-11
o_es_wdb = wilcox(oct, locations = es_wdb) #p = 3.648e-07
o_tm_bg = wilcox(oct, locations = tm_bg) #p = .002717
o_tm_mdb = wilcox(oct, locations = tm_mdb)
o_tm_wda = wilcox(oct, locations = tm_wda)
o_tm_wdb = wilcox(oct, locations = tm_wdb)
o_bg_mdb = wilcox(oct, locations = bg_mdb)
o_bg_wda = wilcox(oct, locations = bg_wda)
o_bg_wdb = wilcox(oct, locations = bg_wdb)
o_mdb_wda = wilcox(oct, locations = mdb_wda)
o_mdb_wdb = wilcox(oct, locations = mdb_wdb)
o_wda_wdb = wilcox(oct, locations = wda_wdb)

ala = clean_all[meta(clean_all)$Species == "alaskensis", ]

a_es_tm = wilcox(ala, locations = es_tm)
a_es_wdb = wilcox(ala, locations = es_wdb)
a_tm_wdb = wilcox(ala, locations = tm_wdb)
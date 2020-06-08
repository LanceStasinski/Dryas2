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
octala = all[!meta(all)$Species == "hybrid"]
meta_octala = meta(octala)
alahyb = all[!meta(all)$Species == "octopetala"]
octhyb = all[!meta(all)$Species == "alaskensis"]

wilcox.test(LMA ~ Species, data = meta_octala, paired = F)#p = 1.5e-13 =dif pops
wilcox.test(LMA ~ Species, data = meta(alahyb), paired =F)#p = .015 = dif pops
wilcox.test(LMA ~ Species, data = meta(octhyb), paired = F)#p = .5977 = NOT DIF

#According to the wilcox test, the LMA for alaskensis is sig. dif. than the LMA 
#of octopetala and hybrid. However, the LMA for hybrids and octopetala are not 
#significantly different. 




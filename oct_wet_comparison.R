################################################################################
#Set up
################################################################################

library("spectrolab")
library("ggplot2")
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

spec1 = readRDS("Clean-up/Clean_spectra/clean_all.rds")
spec = spec1[meta(spec1)$Species == "octopetala"]

################################################################################
#Function
################################################################################

f = function(spec, location){
  spec1 = spec[meta(spec)$Location == location,]
  spec2 = aggregate(spec1, by = meta(spec1)$Name, mean, try_keep_txt(mean))
  spec3 = spec2[,1940]
  specd = data.frame(spec3)
  specd$Location <- location
  colnames(specd) <- c("Reflectance", "Location")
  return(specd) 
}

################################################################################
#By site
################################################################################


es = f(spec, location = "Eagle Summit")
wda = f(spec, location = "Wickersham Dome A")
wdb = f(spec, location = "Wickersham Dome B")

all.df_w_o = Reduce(rbind, list(es, wda, wdb))

################################################################################
#Box Plot
################################################################################

p <- ggplot(all.df_w_o, aes(x=Location, y=Reflectance)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  labs(title="Reflectance at 1940nm - Octopetala(wet)",
       x="Location", y = "Reflectance")+
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("Eagle Summit", "Wickersham Dome A", "Wickersham Dome B"))
theme_classic()

p  

################################################################################
#Tests for normality
################################################################################

shapiro.test(tm$Reflectance) #normal
shapiro.test(es$Reflectance) #not normal
shapiro.test(bg$Reflectance) #normal
shapiro.test(mdb$Reflectance) #normal
shapiro.test(wda$Reflectance) #normal
shapiro.test(wdb$Reflectance) #normal
shapiro.test(all.df$Reflectance) #normal

oct.aov = aov(Reflectance ~ Location, data = all.df_w_o)
summary.aov(oct.aov)

pairwise.t.test(all.df_w_o$Reflectance, all.df_w_o$Location, p.adjust.method = "BH")


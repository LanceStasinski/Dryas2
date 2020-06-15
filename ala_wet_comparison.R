################################################################################
#Set up
################################################################################

library("spectrolab")
library("ggplot2")
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

spec1 = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec = spec1[meta(spec1)$Species == "alaskensis",]

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
wdb = f(spec, location = "Wickersham Dome B")

all.df_a_w = Reduce(rbind, list(es, wdb))

################################################################################
#Box Plot
################################################################################

p <- ggplot(all.df_a_w, aes(x=Location, y=Reflectance)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  labs(title="Vector Normalized Reflectance at 1940nm - Alaskensis(wet)",
       x="Location", y = "Vector Normalized Reflectance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("Eagle Summit", "Wickersham Dome B"))

theme_classic()

p
################################################################################
#Tests for normality
################################################################################


oct.aov = aov(Reflectance ~ Location, data = all.df)
summary.aov(oct.aov)

pairwise.t.test(all.df_a_w$Reflectance, all.df_a_w$Location, p.adjust.method = "BH")

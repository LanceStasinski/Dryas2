library(ggplot2)
library(dplyr)
library(car)
library(vegan)

setwd("~/Dryas Research/Data")

dryas <- read.csv("LMA2.csv")
es = subset(dryas, Location == 2)
tm = subset(dryas, Location == 5)
wdb = subset(dryas, Location == 7)
lma3 = rbind(es, tm, wdb)

dryas3 = lma3

oct = subset(dryas3, Species == 3)
ala = subset(dryas3, Species == 1)
oct.ala = rbind(oct, ala)

################################################################################
#Test for normality of data
################################################################################

shapiro.test(dryas$LMA)

dryas.m = as.matrix(dryas)


mrpp1 = mrpp(dryas.m[,"LMA"], dryas.m[,"Location"])

adonis(LMA ~ Species, data = dryas)
adonis(LMA ~ Location, data = dryas)
adonis(LMA ~ LMA + Location, data = dryas)
adonis(LMA ~ Species * Location, data = dryas)

summary(mrpp1)

################################################################################
#ANOVA all dryas
################################################################################

aov.LMA = aov(LMA ~ Species, data = dryas)
summary.aov(aov.LMA)


TukeyHSD(aov.LMA)

pairwise.t.test(dryas$LMA, dryas$Species, p.adjust.method = "BH")


sp_loc = aov(LMA ~ Species * Location, data = dryas)
summary(sp_loc)

plot(aov.LMA, 1)
plot(aov.LMA, 5)
plot(aov.LMA, 2)

leveneTest(LMA ~ Species, data = dryas)


################################################################################
#Plots - All Sites
################################################################################

ggplot(dryas,aes(x=Species, y=LMA, fill=factor(Species)))+
  geom_boxplot()+ylab("LMA") + 
  stat_summary(fun.y=mean,shape=15,col='blue', size= 2.5,geom='point') +
  ggtitle("LMA per Species - all locations")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))



################################################################################
#ANOVA - 3 sites
################################################################################

aov.LMA3 = aov(LMA ~ Species, data = dryas3)
summary.aov(aov.LMA3)

TukeyHSD(aov.LMA3)
pairwise.t.test(dryas3$LMA, dryas3$Species, p.adjust.method = "BH")

aov.LMA3.loc = aov(LMA ~ Species * Location, data = dryas3)
summary.aov(aov.LMA3.loc)

leveneTest(LMA ~ Species, data = dryas3)

adonis(LMA ~ Species, data = dryas3)

################################################################################
#ANOVA - 3 sites, only oct and ala
################################################################################

aov.octala = aov(LMA ~ Species, data = oct.ala)
summary.aov(aov.octala)

pairwise.t.test(oct.ala$LMA, oct.ala$Species, p.adjust.method = "BH")

aov.octala_loc = aov(LMA ~ Species * Location, data = oct.ala)
summary.aov(aov.octala_loc)

################################################################################
#Plot - 3 sites
################################################################################

ggplot(dryas3,aes(x=Species, y=LMA, fill=factor(Species)))+
  geom_boxplot()+ylab("LMA") + 
  stat_summary(fun.y=mean,shape=15,col='blue', size= 2.5,geom='point') +
  ggtitle("LMA - Cooccurring Sites")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))

################################################################################
#ANOVA - Eagle Summit
################################################################################

eagle.summit<-subset(dryas, Location== "Eagle Summit")
summary(eagle.summit)

aov.ES = aov(LMA ~ Species, data = eagle.summit)
summary.aov(aov.ES)

TukeyHSD(aov.ES)
pairwise.t.test(eagle.summit$LMA, eagle.summit$Species, p.adjust.method = "BH")

leveneTest(LMA ~ Species, data = eagle.summit)

################################################################################
#Plot - Eagle Summit
################################################################################

ggplot(eagle.summit,aes(x=Species, y=LMA, fill=factor(Species)))+
  geom_boxplot()+ylab("LMA") + 
  stat_summary(fun.y=mean,shape=15,col='blue', size= 2.5,geom='point') +
  ggtitle("Eagle Summit")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))

################################################################################
#ANOVA - Twelve Mile
################################################################################

twelve.mile<-subset(dryas, Location== "Twelve Mile")

aov.TM = aov(LMA ~ Species, data = twelve.mile)
summary.aov(aov.TM)

TukeyHSD(aov.TM)
pairwise.t.test(twelve.mile$LMA, twelve.mile$Species, p.adjust.method = "BH")

leveneTest(LMA ~ Species, data = twelve.mile)


################################################################################
#Plot - Twelve Mile
################################################################################

twelve.mile<-subset(dryas, Location== "Twelve Mile")

ggplot(twelve.mile,aes(x=Species, y=LMA, fill=factor(Species)))+
  geom_boxplot()+ylab("LMA") + 
  stat_summary(fun.y=mean,shape=15,col='blue', size= 2.5,geom='point') +
  ggtitle("Twelve Mile")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))

################################################################################
#ANOVA - Wickersham Dome B
################################################################################

wickersham.b<-subset(dryas, Location== "Wickersham Dome B")

aov.WDB = aov(LMA ~ Species, data = wickersham.b)
summary.aov(aov.WDB)

TukeyHSD(aov.WDB)
pairwise.t.test(wickersham.b$LMA, wickersham.b$Species, p.adjust.method = "BH")

leveneTest(LMA ~ Species, data = wickersham.b)

################################################################################
#Plot - Wickersham Dome B
################################################################################

wickersham.b<-subset(dryas, Location== "Wickersham Dome B")

ggplot(wickersham.b,aes(x=Species, y=LMA, fill=factor(Species)))+
  geom_boxplot()+ylab("LMA") + 
  stat_summary(fun.y=mean,shape=15,col='blue', size= 2.5,geom='point') +
  ggtitle("Wickersham Dome B")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=10))


################################################################################
#D. octopetala
################################################################################

oct = subset(dryas, Species == "octopetala")

oct_loc = aov(LMA ~ Location, data = oct)
summary(oct_loc)

pairwise.t.test(oct$LMA, oct$Location, p.adjust.method = "BH")

################################################################################
#D. alaskensis
################################################################################

ala = subset(dryas, Species == "alaskensis")

ala_loc = aov(LMA ~ Location, data = ala)
summary(ala_loc)

pairwise.t.test(ala$LMA, ala$Location, p.adjust.method = "BH")

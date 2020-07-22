#PROSPECT

library(spectrolab)
library(ggplot2)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

spec = readRDS("Clean-up/Vector_normalized/all_vn.rds")
log = read.csv("log.csv")

meta(spec) <- log
spec = spec[!meta(spec)$Species_ID == "DX",]
data = meta(spec)

par(mar = c(4,4,2,1))
plot(data$N.mu, main = "Number of leaf layers (N) predicted by PROSPECT", xlab = "Sample number", 
     ylab = "N")


oct = subset(data, data$Species_ID == "DO")
ala = subset(data, data$Species_ID == "DA")

p <- ggplot(data, aes(x=Species_ID, y=Cw.mu)) + 
  geom_boxplot(fill='#A4A4A4', color="black")+
  labs(title="Cw.mu per species",
       x="Species", y = "Cw.mu") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits=c("DA", "DO", "DX")) 


p

boxplot(data$LMA, data$Cm.mu, main = "LMA vs Cm.mu", names = c("LMA", "Cm.mu"))

boxplot(data$LMA, data$Cm.med, main = "LMA vs Cm.med", names = c("LMA", "Cm.med"))

boxplot(data$LMA, data$Cm.q975, main = "LMA vs Cm.q975", names = c("LMA", "Cm.q975"))

lma.mu = mean(data$LMA)
t.test(data$Cm.q975, mu = lma.mu, alternative = "two.sided")





plot(sort(data$Cm.mu, decreasing = F), type = "l", col = "red", lwd = 1, ylim = c(0, 0.025))
lines(sort(data$LMA, decreasing = F), col="blue",lty=1, lwd = 1)

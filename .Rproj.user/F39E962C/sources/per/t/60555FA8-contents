#Morphological trait analysis and plotting
################################################################################
#Set up
################################################################################

library(corrplot)
library(plyr)
library(vegan)
library(ggplot2)
library(ggpubr)
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

leaf = read.csv('morphology_2.csv', stringsAsFactors = F)
values = c("parent",'parent','hybrid')
leaf$Species = as.factor(leaf$Species)
leaf$taxa <- values[leaf$Species]

################################################################################
#Manova
################################################################################
#parents vs hybrids
gm = leaf$Glandular.Midvien
rs = leaf$Rusty.Scales
l = leaf$Length
at = leaf$Adaxial.tomentum
traits = Reduce(cbind, list(gm,rs,l,at))

m = manova(traits ~ Species, data = leaf)
summary(m)

#DA vs DO
nohyb = leaf[!leaf$Species == 'DX',]
gm = nohyb$Glandular.Midvien
rs = nohyb$Rusty.Scales
l = nohyb$Length
at = nohyb$Adaxial.tomentum
traits = Reduce(cbind, list(gm,rs,l,at))

summary(manova(traits ~ Species, data = nohyb))

#Hybrid vs DA
nodo = leaf[!leaf$Species == 'DO',]
gm = nodo$Glandular.Midvien
rs = nodo$Rusty.Scales
l = nodo$Length
at = nodo$Adaxial.tomentum
traits = Reduce(cbind, list(gm,rs,l,at))

summary(manova(traits ~ Species, data = nodo))

#hybrid vs DO
noda = leaf[!leaf$Species == 'DA',]
gm = noda$Glandular.Midvien
rs = noda$Rusty.Scales
l = noda$Length
at = noda$Adaxial.tomentum
traits = Reduce(cbind, list(gm,rs,l,at))

summary(manova(traits ~ Species, data = noda))

################################################################################
#Mantel test
################################################################################
traits = leaf[, -c(1,2,3,9,10,11)]
parent = traits[!traits$Species == 'DX',]
parent = parent[,-5]
hyb = traits[traits$Species == 'DX',]
hyb = hyb[,-5]

cor.p = cor(parent, method = 'pearson')
cor.h = cor(hyb, method = 'pearson')

mantel.t = mantel(cor.p, cor.h, method = 'pearson', permutations = 9999)
mantel.t
################################################################################
#Trait correlations
################################################################################
#data sets
traits = leaf[, -c(1,2,3,9,10,11)]
parent = traits[!traits$Species == 'DX',]
parent = parent[,-5]
hyb = traits[traits$Species == 'DX',]
hyb = hyb[,-5]


#parent trait correlations
p.c = cor(parent, method = 'pearson')
colnames(p.c) <- c('Glands', 'Scales', 'Length', 'Tomentum')
rownames(p.c) <- c('Glands', 'Scales', 'Length', 'Tomentum')


#hybrid trait correlations
hyb.c = cor(hyb, method = 'pearson')
colnames(hyb.c) <-c('Glands', 'Scales', 'Length', 'Tomentum')
rownames(hyb.c) <- c('Glands', 'Scales', 'Length', 'Tomentum')

################################################################################
#Plot traits by species
################################################################################
#length
bp = ggplot(leaf, aes(x = Species, y = Length))+ ylab('Length (mm)')+
        geom_boxplot() +
        ggtitle('Leaf Length')+
        scale_x_discrete(labels=c('DAK', 'DAJ', 'DX'))+
        theme(plot.title = element_text(hjust = .5))

#Midvein total
midvein = ggplot(leaf,
               aes(x = factor(Species,
                              levels = c('DA', 'DO', 'DX'),
                              labels = c('DAK', 'DAJ', 'DX')),
                   fill = factor(Midvein,
                                 levels = c(0,1,2,3),
                                 labels = c('Neither', 'Glands', 'Scales',
                                            'Glands + Scales')))) +
  geom_bar(position = 'fill') +
  scale_y_continuous(breaks = seq(0,1,.2)) +
  scale_fill_grey(start = .9, end = .1) +
  labs(y = "Proportion",
       fill = 'Glands or Scales',
       x = 'Species') +
  theme_minimal()+
  ggtitle('Abaxial Midvein Morphology')+
  theme(plot.title = element_text(hjust = .5))

#Midvein Glands
gland = ggplot(leaf,
       aes(x = factor(Species,
                      levels = c('DA', 'DO', 'DX')),
           fill = factor(Glandular.Midvien,
                         levels = c(0,1),
                         labels = c('Absent', 'Present')))) +
  geom_bar(position = 'fill') +
  scale_y_continuous(breaks = seq(0,1,.2)) +
  scale_fill_grey(start = .8, end = .2) +
  labs(y = "Proportion",
       fill = 'Glands',
       x = 'Species') +
  theme_minimal()+
  ggtitle('Midvein Glandular Trichomes')+
  theme(plot.title = element_text(hjust = .5))

#Midvein Scales
scale = ggplot(leaf,
       aes(x = factor(Species,
                      levels = c('DAK', 'DAJ', 'DX')),
           fill = factor(Rusty.Scales,
                         levels = c(0,1),
                         labels = c('Absent', 'Present')))) +
        geom_bar(position = 'fill') +
        scale_y_continuous(breaks = seq(0,1,.2)) +
        scale_fill_grey(start = .8, end = .2) +
        labs(y = "Proportion",
             fill = 'Scales',
             x = 'Species') +
        theme_minimal()+
        ggtitle('Midvein Scales')+
        theme(plot.title = element_text(hjust = .5))

#has both scales and glands
scale_gland = ggplot(leaf,
               aes(x = factor(Species,
                              levels = c('DAK', 'DAJ', 'DX')),
                   fill = factor(Scales.and.glands,
                                 levels = c(0,1),
                                 labels = c('Not concurrent', 'Concurrent')))) +
  geom_bar(position = 'fill') +
  scale_y_continuous(breaks = seq(0,1,.2)) +
  scale_fill_grey(start = .8, end = .2) +
  labs(y = "Proportion",
       fill = 'Scales and Glands',
       x = 'Species') +
  theme_minimal()+
  ggtitle('Pressence of scales and glands')+
  theme(plot.title = element_text(hjust = .5))

#tomentum
tom = ggplot(leaf,
       aes(x = factor(Species,
                      levels = c('DA', 'DO', 'DX'),
                      labels = c('DAK', 'DAJ', 'DX')),
           fill = factor(Adaxial.tomentum,
                         levels = c(1,2,3),
                         labels = c('Sparse', 'Moderate', 'Dense')))) +
        geom_bar(position = 'fill') +
        scale_y_continuous(breaks = seq(0,1,.2)) +
        scale_fill_grey(start = .8, end = .2) +
        labs(y = "Proportion",
             fill = 'Tomentum',
             x = 'Species') +
        theme_minimal()+
        ggtitle('Adaxial Tomentum')+
        theme(plot.title = element_text(hjust = .5))

ggarrange(bp, midvein, tom, ncol = 2, nrow = 2)

#corrplots
par(mfrow = c(1,2))
corrplot(p.c, type = "lower", diag = F, cl.cex = 1, addCoef.col = 'black',
         tl.col = 'black', cl.length = 5)
mtext('Parent', side = 2, line = -1, at = 2, cex = 1.25)
corrplot(hyb.c, type = 'lower', diag = F, cl.cex = 1, addCoef.col = 'black',
         tl.col = 'black', cl.length = 5)
mtext('Hybrid', side = 2, line = -1, at = 2, cex = 1.25)

################################################################################
#Bootstrapping - a different way to compare traits - not used in the paper
################################################################################
library(boot)
#scale vs gland
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,1],d[,2], method = 'spearman'))
}
p.scale_gland = boot(parent, statistic = correlation, R = 9999)
p.sg_ci = boot.ci(p.scale_gland, conf = .95)
h.scale_gland = boot(hyb, statistic = correlation, R = 9999)
h.sg_ci = boot.ci(h.scale_gland, conf = .95)
n.scale_gland = boot(null.traits, statistic = correlation, R = 9999)
n.sg_ci = boot.ci(n.scale_gland, conf = .95)

#gland vs length
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,1],d[,3], method = "spearman"))
}
p.gl = boot(parent, statistic = correlation, R =9999)
p.gl_ci = boot.ci(p.gl, conf = .95)
h.gl = boot(hyb, statistic = correlation, R =9999)
h.gl_ci = boot.ci(h.gl, conf = .95)
n.gl = boot(null.traits, statistic = correlation, R = 9999)
n.gl_ci = boot.ci(n.gl, conf = .95)

#gland vs tomentum
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,1],d[,4],method = 'spearman'))
}
p.gt = boot(parent, statistic = correlation, R =9999)
p.gt_ci = boot.ci(p.gt, conf = .95)
h.gt = boot(hyb, statistic = correlation, R =9999)
h.gt_ci = boot.ci(h.gt, conf = .95)
n.gt = boot(null.traits, statistic = correlation, R = 9999)
n.gt_ci = boot.ci(n.gt, conf = .95)
#scale vs length
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,2],d[,3], method = "spearman"))
}
p.sl = boot(parent, statistic = correlation, R =9999)
p.sl_ci = boot.ci(p.sl, conf = .95)
h.sl = boot(hyb, statistic = correlation, R =9999)
h.sl_ci = boot.ci(h.sl, conf = .95)
n.sl = boot(null.traits, statistic = correlation, R = 9999)
n.sl_ci = boot.ci(n.sl, conf = .95)

#scale vs tomentum
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,2],d[,4],method = 'spearman'))
}
p.st = boot(parent, statistic = correlation, R =9999)
p.st_ci = boot.ci(p.st, conf = .95)
h.st = boot(hyb, statistic = correlation, R =9999)
h.st_ci = boot.ci(h.st, conf = .95)
n.st = boot(null.traits, statistic = correlation, R = 9999)
n.st_ci = boot.ci(n.st, conf = .95)

#length vs tomentum
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,3],d[,4], method = "spearman"))
}
p.lt = boot(parent, statistic = correlation, R =9999)
p.lt_ci = boot.ci(p.lt, conf = .95)
h.lt = boot(hyb, statistic = correlation, R =9999)
h.lt_ci = boot.ci(h.lt, conf = .95)
n.lt = boot(null.traits, statistic = correlation, R = 9999)
n.lt_ci = boot.ci(n.lt, conf = .95)


################################################################################
#Plot bootstrap distributions separately
################################################################################

par(mfrow = c(2,6))
#Parents
hist(p.scale_gland$t, main = "Scales vs Glands", xlab = NA, 
     ylab = "Frequency - Parent Species", xlim = c(-1, 0))
abline(v = p.sg_ci$bca[4], lty = 2)
abline(v = p.sg_ci$bca[5], lty = 2)
abline(v= p.scale_gland$t0, lty = 1, lwd = 2, col = "red")

hist(p.gl$t, main = 'Glands vs Length', xlab = NA, ylab = NA,
     xlim = c(-.2, 1))
abline(v = p.gl_ci$bca[5], lty =2)
abline(v = p.gl_ci$bca[4], lty =2)
abline(v= p.gl$t0, lty = 1, lwd = 2, col = "red")

hist(p.gt$t, main = 'Glands vs Tomentum', xlab = NA, ylab = NA,
     xlim = c(-1,.2))
abline(v = p.gt_ci$bca[5], lty =2)
abline(v = p.gt_ci$bca[4], lty =2)
abline(v= p.gt$t0, lty = 1, lwd = 2, col = "red")

hist(p.sl$t, main = 'Scales vs Length', xlab = NA, ylab = NA,
     xlim = c(-1,.4))
abline(v = p.sl_ci$bca[5], lty =2)
abline(v = p.sl_ci$bca[4], lty =2)
abline(v= p.sl$t0, lty = 1, lwd = 2, col = "red")

hist(p.st$t, main = 'Scales vs Tomentum', xlab = NA, ylab = NA,
     xlim = c(-1,1))
abline(v = p.st_ci$bca[5], lty =2)
abline(v = p.st_ci$bca[4], lty =2)
abline(v= p.st$t0, lty = 1, lwd = 2, col = "red")

hist(p.lt$t, main = 'Length vs Tomentum', xlab = NA, ylab = NA,
     xlim = c(-1,1))
abline(v = p.lt_ci$bca[5], lty =2)
abline(v = p.lt_ci$bca[4], lty =2)
abline(v= p.lt$t0, lty = 1, lwd = 2, col = "red")

#hybrids
hist(h.scale_gland$t, main = NA, xlab = NA, 
     ylab = "Frequency - Hybrids", xlim = c(-1, 0))
abline(v = h.sg_ci$bca[4], lty = 2)
abline(v = h.sg_ci$bca[5], lty = 2)
abline(v=h.scale_gland$t0, lty = 1, lwd = 2, col = "red")

hist(h.gl$t, main = NA, xlab = NA, ylab = NA, xlim = c(-.2,1))
abline(v = h.gl_ci$bca[5], lty =2)
abline(v = h.gl_ci$bca[4], lty =2)
abline(v=h.gl$t0, lty = 1, lwd = 2, col = "red")

hist(h.gt$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,.2))
abline(v = h.gt_ci$bca[5], lty =2)
abline(v = h.gt_ci$bca[4], lty =2)
abline(v=h.gt$t0, lty = 1, lwd = 2, col = "red")

hist(h.sl$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,.4))
abline(v = h.sl_ci$bca[5], lty =2)
abline(v = h.sl_ci$bca[4], lty =2)
abline(v=h.sl$t0, lty = 1, lwd = 2, col = "red")

hist(h.st$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,1))
abline(v = h.st_ci$bca[5], lty =2)
abline(v = h.st_ci$bca[4], lty =2)
abline(v=h.st$t0, lty = 1, lwd = 2, col = "red")

hist(h.lt$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,1))
abline(v = h.lt_ci$bca[5], lty =2)
abline(v = h.lt_ci$bca[4], lty =2)
abline(v= h.lt$t0, lty = 1, lwd = 2, col = "red")

################################################################################
#Plot bootstrap distributions side by side
################################################################################

par(mfrow = c(5,2))
boxplot(leaf$Length~leaf$Species, ylab = 'Species', xlab = 'Leaf Length (mm)',
        main = 'Leaf length', col = c("#00B0F6","#F8766D",'gray30'), 
        horizontal = T, las = 1)
barplot(glands, main="Glandular midvein",
        xlab= NA, names = c("Absent", "Present"),
        ylab = "Counts", col=c("#00B0F6","#F8766D",'gray30'), beside=TRUE)
barplot(scales, main = "Midvein scales", 
        xlab = NA, names = c("Absent", "Present"),
        ylab = "Counts", col = c("#00B0F6", "#F8766D", "gray30"), beside=TRUE)
barplot(adaxial, main = "Adaxial tomentum", 
        xlab = NA, names = c('Sparse', 'Moderate', 'Dense'),
        ylab = "Counts", col = c("#00B0F6", "#F8766D", "gray30"), beside=TRUE)

cp <- rgb(210,0,255, max = 255, alpha = 70)
cp2 = rgb(210,0,255, max = 255, alpha = 255)
ch <- rgb(0, 0, 0, max = 255, alpha = 70)

#Scale-gland
hist(p.scale_gland$t, main = "Scales vs Glands", xlab = 'Rho', 
     ylab = 'Frequency', xlim = c(-1, 0), col = cp)
abline(v = p.sg_ci$bca[4], lty = 2, col = cp2, lwd = 2)
abline(v = p.sg_ci$bca[5], lty = 2, col = cp2, lwd = 2)
abline(v= p.scale_gland$t0, lty = 1, lwd = 3, col = cp2)

hist(h.scale_gland$t, main = NA, xlab = NA, 
     ylab = NA, xlim = c(-1, 0), add = T, col = ch)
abline(v = h.sg_ci$bca[4], lty = 2, col = "black", lwd = 2)
abline(v = h.sg_ci$bca[5], lty = 2, col = "black", lwd = 2)
abline(v=h.scale_gland$t0, lty = 1, lwd = 3, col = "black")


#Gland-length
hist(p.gl$t, main = 'Glands vs Length', xlab = 'Rho', ylab = 'Frequency',
     xlim = c(-.2, 1), col = cp)
abline(v = p.gl_ci$bca[5], lty =2, col = cp2, lwd = 2)
abline(v = p.gl_ci$bca[4], lty =2, col = cp2, lwd = 2)
abline(v= p.gl$t0, lty = 1, lwd = 3, col = cp2)

hist(h.gl$t, main = NA, xlab = NA, ylab = NA, xlim = c(-.2,1), add = T, col = ch)
abline(v = h.gl_ci$bca[5], lty =2, col = "black", lwd = 2)
abline(v = h.gl_ci$bca[4], lty =2, col = "black", lwd = 2)
abline(v=h.gl$t0, lty = 1, lwd = 3, col = "black")

#Gland-tomentum
hist(p.gt$t, main = 'Glands vs Tomentum', xlab = 'Rho', ylab = 'Frequency',
     xlim = c(-1,.2), col = cp)
abline(v = p.gt_ci$bca[5], lty =2, col = cp2, lwd = 2)
abline(v = p.gt_ci$bca[4], lty =2, col = cp2, lwd = 2)
abline(v= p.gt$t0, lty = 1, lwd = 3, col = cp2)

hist(h.gt$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,.2), add = T, col = ch)
abline(v = h.gt_ci$bca[5], lty =2, col = "black", lwd = 2)
abline(v = h.gt_ci$bca[4], lty =2, col = "black", lwd = 2)
abline(v=h.gt$t0, lty = 1, lwd = 3, col = "black")

#Scales-length
hist(p.sl$t, main = 'Scales vs Length', xlab = 'Rho', ylab = 'Frequency',
     xlim = c(-1,.4), col = cp)
abline(v = p.sl_ci$bca[5], lty =2, col = cp2, lwd = 2)
abline(v = p.sl_ci$bca[4], lty =2, col = cp2, lwd = 2)
abline(v= p.sl$t0, lty = 1, lwd = 3, col = cp2)

hist(h.sl$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,.4), add = T, col = ch)
abline(v = h.sl_ci$bca[5], lty =2, col = "black", lwd = 2)
abline(v = h.sl_ci$bca[4], lty =2, col = "black", lwd = 2)
abline(v=h.sl$t0, lty = 1, lwd = 3, col = "black")

#Scales-tomentum
hist(p.st$t, main = 'Scales vs Tomentum', xlab = 'Rho', ylab = 'Frequency',
     xlim = c(-.5,1), col = cp)
abline(v = p.st_ci$bca[5], lty =2, col = cp2, lwd = 2)
abline(v = p.st_ci$bca[4], lty =2, col = cp2, lwd = 2)
abline(v= p.st$t0, lty = 1, lwd = 3, col = cp2)

hist(h.st$t, main = NA, xlab = NA, ylab = NA, xlim = c(-.5,1), add = T, col = ch)
abline(v = h.st_ci$bca[5], lty =2, col = "black", lwd = 2)
abline(v = h.st_ci$bca[4], lty =2, col = "black", lwd = 2)
abline(v=h.st$t0, lty = 1, lwd = 3, col = "black")

#Length-tomentum
hist(p.lt$t, main = 'Length vs Tomentum', xlab = 'Rho', ylab = 'Frequency',
     xlim = c(-1,1), col = cp)
abline(v = p.lt_ci$bca[5], lty =2, col = cp2, lwd = 2)
abline(v = p.lt_ci$bca[4], lty =2, col = cp2, lwd = 2)
abline(v= p.lt$t0, lty = 1, lwd = 3, col = cp2)

hist(h.lt$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,1), add = T, col = ch)
abline(v = h.lt_ci$bca[5], lty =2, col = "black", lwd = 2)
abline(v = h.lt_ci$bca[4], lty =2, col = "black", lwd = 2)
abline(v= h.lt$t0, lty = 1, lwd = 3, col = "black")








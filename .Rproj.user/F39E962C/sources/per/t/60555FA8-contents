library(polycor)
library(corrplot)
library(plyr)
#Morphology analysis

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

leaf = read.csv('morphology.csv', stringsAsFactors = F)

################################################################################
#Plot traits by species
################################################################################

glands = table(leaf$Species, leaf$Glandular.Midvien)
scales = table(leaf$Species, leaf$Rusty.Scales)
adaxial = table(leaf$Species, leaf$Adaxial.tomentum)

par(mfrow = c(2,2))
boxplot(leaf$Length~leaf$Species, xlab = NA, ylab = 'Leaf Length',
        main = 'Leaf length')
barplot(glands, main="Glandular midvein",
        xlab= NA, names = c("Absent", "Present"),
        ylab = "Counts", col=c("#00B0F6","#F8766D",'gray30'), beside=TRUE)
barplot(scales, main = "Midvein scales", 
        xlab = NA, names = c("Absent", "Present"),
        ylab = "Counts", col = c("#00B0F6", "#F8766D", "gray30"), beside=TRUE)

barplot(adaxial, main = "Adaxial tomentum", 
        xlab = NA, names = c('Sparse', 'Moderate', 'Dense'),
        ylab = "Counts", col = c("#00B0F6", "#F8766D", "gray30"), beside=TRUE)

################################################################################
#corrplots
################################################################################
traits = leaf[, -c(1,2,3,9,10)]

#parent trait correlations
parent = traits[!traits$Species == 'DX',]
parent = parent[,-5]
p.c = cor(parent, method = 'spearman')
corrplot(p.c, type = "lower", diag = F)

#hybrid trait correlations
hyb = traits[traits$Species == 'DX',]
hyb = hyb[,-5]
hyb.c = cor(hyb, method = 'spearman')
corrplot(hyb.c, type = 'lower', diag = F)

################################################################################
#Bootstrapping
################################################################################
library(boot)
#scale vs gland
correlation <- function(data, indices) {
        d <- data[indices,]
        return(polychor(d[,1],d[,2]))
}
p.scale_gland = boot(parent, statistic = correlation, R = 1000)
p.sg_ci = boot.ci(p.scale_gland, conf = .95)
h.scale_gland = boot(hyb, statistic = correlation, R = 1000)
h.sg_ci = boot.ci(h.scale_gland, conf = .95)

#gland vs length
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,1],d[,3], method = "spearman"))
}
p.gl = boot(parent, statistic = correlation, R =1000)
p.gl_ci = boot.ci(p.gl, conf = .95)
h.gl = boot(hyb, statistic = correlation, R =1000)
h.gl_ci = boot.ci(h.gl, conf = .95)

#gland vs tomentum
correlation <- function(data, indices) {
        d <- data[indices,]
        return(polychor(d[,1],d[,4]))
}
p.gt = boot(parent, statistic = correlation, R =1000)
p.gt_ci = boot.ci(p.gt, conf = .95)
h.gt = boot(hyb, statistic = correlation, R =1000)
h.gt_ci = boot.ci(h.gt, conf = .95)


#scale vs length
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,2],d[,3], method = "spearman"))
}
p.sl = boot(parent, statistic = correlation, R =1000)
p.sl_ci = boot.ci(p.sl, conf = .95)
h.sl = boot(hyb, statistic = correlation, R =1000)
h.sl_ci = boot.ci(h.sl, conf = .95)

#scale vs tomentum
correlation <- function(data, indices) {
        d <- data[indices,]
        return(polychor(d[,2],d[,4]))
}
p.st = boot(parent, statistic = correlation, R =1000)
p.st_ci = boot.ci(p.st, conf = .95)
h.st = boot(hyb, statistic = correlation, R =1000)
h.st_ci = boot.ci(h.st, conf = .95)

#length vs tomentum
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,3],d[,4], method = "spearman"))
}
p.lt = boot(parent, statistic = correlation, R =1000)
p.lt_ci = boot.ci(p.lt, conf = .95)
h.lt = boot(hyb, statistic = correlation, R =1000)
h.lt_ci = boot.ci(h.lt, conf = .95)

#Plot

par(mfrow = c(2,6))
#Parents
hist(p.scale_gland$t, main = "Scales vs Glands", xlab = NA, 
     ylab = "Frequency - Parent Species", xlim = c(-1, 0))
abline(v = p.sg_ci$bca[4], lty = 2)
abline(v = p.sg_ci$bca[5], lty = 2)


hist(p.gl$t, main = 'Glands vs Length', xlab = NA, ylab = NA,
     xlim = c(-.2, 1))
abline(v = p.gl_ci$bca[5], lty =2)
abline(v = p.gl_ci$bca[4], lty =2)

hist(p.gt$t, main = 'Glands vs Tomentum', xlab = NA, ylab = NA,
     xlim = c(-1,.2))
abline(v = p.gt_ci$bca[5], lty =2)
abline(v = p.gt_ci$bca[4], lty =2)

hist(p.sl$t, main = 'Scales vs Length', xlab = NA, ylab = NA,
     xlim = c(-1,.4))
abline(v = p.sl_ci$bca[5], lty =2)
abline(v = p.sl_ci$bca[4], lty =2)

hist(p.st$t, main = 'Scales vs Tomentum', xlab = NA, ylab = NA,
     xlim = c(-1,1))
abline(v = p.st_ci$bca[5], lty =2)
abline(v = p.st_ci$bca[4], lty =2)

hist(p.lt$t, main = 'Length vs Tomentum', xlab = NA, ylab = NA,
     xlim = c(-1,1))
abline(v = p.lt_ci$bca[5], lty =2)
abline(v = p.lt_ci$bca[4], lty =2)

#hybrids
hist(h.scale_gland$t, main = NA, xlab = NA, 
     ylab = "Frequency - Hybrids", xlim = c(-1, 0))
abline(v = h.sg_ci$bca[4], lty = 2)
abline(v = h.sg_ci$bca[5], lty = 2)
abline(v=mean(h.scale_gland$t), lty = 1, lwd = 2, col = "red")

hist(h.gl$t, main = NA, xlab = NA, ylab = NA, xlim = c(-.2,1))
abline(v = h.gl_ci$bca[5], lty =2)
abline(v = h.gl_ci$bca[4], lty =2)
abline(v=mean(h.gl$t), lty = 1, lwd = 2, col = "red")

hist(h.gt$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,.2))
abline(v = h.gt_ci$bca[5], lty =2)
abline(v = h.gt_ci$bca[4], lty =2)
abline(v=mean(h.gt$t), lty = 1, lwd = 2, col = "red")

hist(h.sl$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,.4))
abline(v = h.sl_ci$bca[5], lty =2)
abline(v = h.sl_ci$bca[4], lty =2)
abline(v=mean(h.sl$t), lty = 1, lwd = 2, col = "red")

hist(h.st$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,1))
abline(v = h.st_ci$bca[5], lty =2)
abline(v = h.st_ci$bca[4], lty =2)
abline(v=mean(h.st$t), lty = 1, lwd = 2, col = "red")

hist(h.lt$t, main = NA, xlab = NA, ylab = NA, xlim = c(-1,1))
abline(v = h.lt_ci$bca[5], lty =2)
abline(v = h.lt_ci$bca[4], lty =2)
abline(v=mean(h.lt$t), lty = 1, lwd = 2, col = "red")




#t tests



t.test(p.scale_gland$t, h.scale_gland$t)
t.test(p.gl$t, h.gl$t)
t.test(p.gt$t, h.gt$t)
t.test(p.sl$t, h.sl$t)
t.test(p.st$t, h.st$t)
t.test(p.lt$t, h.lt$t)


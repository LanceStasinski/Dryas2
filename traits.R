#Morphology analysis

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

leaf = read.csv('morphology.csv', stringsAsFactors = F)

#Abaxial Tomentum

count = table(leaf$Species, leaf$Abaxial.tomentum)

barplot(count, main="Abaxial tomentum by species",
        xlab="Abaxial tomentum", ylab = "Counts",
        col=c("Blue","red",'black'),
        legend = rownames(count), beside=TRUE)


#Plot
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

#correlations
parent = leaf[!leaf$Species == 'DX',]
cor.parent = cor(parent$Glandular.Midvien, parent$Rusty.Scales)
hyb = leaf[leaf$Species == 'DX',]
cor.hyb = cor(hyb$Glandular.Midvien, hyb$Rusty.Scales)

library(boot)
correlation <- function(data, indices) {
        d <- data[indices,]
        return(cor(d[,4],d[,5]))
}
bt = boot(parent, statistic = correlation, R = 1000)
boot.ci(bt, conf = .95)


bt2 = boot(hyb, statistic = correlation, R = 1000)
boot.ci(bt2, conf = .95)

par(mfrow = c(2,1))
hist(bt$t, main = "Parent Species", xlab = NA,
     xlim = c(-1, 0))
abline(v = -.966, lty = 2)
abline(v = -.772, lty = 2)
hist(bt2$t, main = "Hybrids", xlim = c(-1, 0),
     xlab = "Correlation between midvein glands and midvein scales")
abline(v = -.7385, lty = 2)
abline(v = -.2708, lty = 2)


t.test(bt$t, bt2$t)

da = leaf[!leaf$DA == 'NA',]
cor(da$Glandular.Midvien, da$DA)



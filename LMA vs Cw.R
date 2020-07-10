library(spectrolab)

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
# Fit PLS_DA model all
################################################################################
#data
data1 = read.csv("log.csv")
spec_all = readRDS("Clean-up/Vector_normalized/all_vn.rds")

meta(spec_all) <- data1

spec_all = spec_all[!meta(spec_all)$Cw.mu == "NA",]
data = meta(spec_all)


plot(data$LMA, data$Cw.mu, ylab = "Cw.mu", xlab = "LMA", main = "LMA vs Cw.mu")
abline(lm(Cw.mu ~ LMA, data = data), col = "blue")

par(mar = c(5, 5, 3, 5))
plot(data$LMA, type ="p", pch = 16, ylab = "LMA", ylim = c(0, 0.025),
     main = "LMA and Cw.mu", xlab = "Sample #",
     col = "blue")
par(new = TRUE)
plot(data$Cw.mu * 4.876836, type = "p", xaxt = "n", ylim = c(0, 0.025), yaxt = "n",
     ylab = "", xlab = "", col = "red", pch = 17)
axis(side = 4)
mtext("Cw.mu", side = 4, line = 3)
legend("topright", c("LMA", "Cw.mu"),
       col = c("blue", "red"), pch = c(16,17))

cw = mean(data$Cw.mu)
lma = mean(data$LMA)



#PLSDA
#RDS files for cleaned data can be downloaded from:
#FM_Polar_Studies > Dryas_Spectral_Analyses > Clean-up 
#or generated via 01_clean_up.R
################################################################################
#Set up
################################################################################

library(mixOmics)
library(spectrolab)
library(plotrix)
library(rgl)
library(grDevices)
################################################################################
#Set working directory to folder containing downloaded rds files
################################################################################

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
# Fit PLS_DA model all
################################################################################
#data
spec_all = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_es = spec_all[meta(spec_all)$Location == "Eagle Summit",]
spec_wdb = spec_all[meta(spec_all)$Location == "Wickersham Dome B",]
spec_tm = spec_all[meta(spec_all)$Location == "Twelve Mile",]
spec_all = Reduce(spectrolab::combine, list(spec_es, spec_wdb, spec_tm))

names(spec_all) = meta(spec_all)$Species_ID
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 50)

perf.plot_species = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
saveRDS(perf.plot_species, "Figures/perf plots/perf.plot_species.rds")
###ncomp = 16
plotIndiv(plsda.fit, title = "", comp = c(1,2,5), legend = TRUE, 
          style = "3d", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 40)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,26] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)

par(mar = c(2, 5, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1, line = -1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)
mtext("Predict Species, 3 sites, 26 comps, 0.2051 BER", side = 3, line = 1.5)

library(stringr)

str_locate_all(pattern ='DX', prediction)
h = c(76,77,78,84,85)
prediction[h]
hspec = c(217,221,223,252,253)


r = resp[test]
str_locate_all(pattern ='DX', r)
a = c(33,34,35,76,77,78,79,83,84,85,132,133,134)
prediction[a]
aspec = c(81,225,251,377)
ospec = c(83,84,373,380)

pred.as.hyb = spec_all[hspec,]
pred.hyb.DA = meta(pred.as.hyb)$DA

pred.as.ala = spec_all[aspec,]
pred.ala.DA = meta(pred.as.ala)$DA

pred.as.oct = spec_all[ospec,]
pred.oct.DA = meta(pred.as.oct)$DA

par(mar = c(4, 4, 3, 1), oma = c(2, 4, 3, 2))
plot(y = c(1:5), x = pred.hyb.DA, type = 'p', pch = 16, col = 'black', cex = 1.2,
     ylab = "Sample", xlab = "Proportion of DA ancestry", xlim = c(.3,.7),
     main = "Ancestry from Alaskensis - samples predicted as hybrids")
points(y = c(1:4), x = pred.oct.DA, type = 'p', pch = 15, cex = 1.5, col = "blue")
points(y = c(1:4), x = pred.ala.DA, type = 'p', pch = 17, cex = 1.0, col = "orange")
abline(v = .5, lty = 2, lwd = 1)
legend("topright", legend = c("Predicted DX", "Predicted DO", "Predicted DA"),
       pch = c(16,15,17), col = c("black", "blue", "orange"))






################################################################################
# Fit PLS_DA model all
################################################################################


spec_all = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_all = spec_all[meta(spec_all)$Species_ID == "DX",]
names(spec_all) = meta(spec_all)$Location
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 20)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 50)

perf.plot_species = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
                         legend.position = "horizontal")
saveRDS(perf.plot_species, "Figures/perf plots/perf.plot_species.rds")
###ncomp = 16
plotIndiv(plsda.fit, title = "", comp = c(1,2,5), legend = TRUE, 
          style = "3d", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 10)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,7] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)


#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)

################################################################################
# PLSDA predict species and location
################################################################################

#data
spec_all = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec_all = spec_all[!meta(spec_all)$sp_loc == "NaN",]
spec_all = spec_all[!meta(spec_all)$Species_ID == "DX",]
names(spec_all) = meta(spec_all)$sp_loc
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 50)

perf.plot_species = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
                         legend.position = "horizontal")
saveRDS(perf.plot_species, "Figures/perf plots/perf.plot_species.rds")
###ncomp = 16
plotIndiv(plsda.fit, title = "", comp = c(1,2,5), legend = TRUE, 
          style = "3d", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 30)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,24] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)


#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = .75)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)


################################################################################
# PLSDA simulate airborne sensor
################################################################################

#data
spec_all = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec_all = spec_all[meta(spec_all)$Location == "Eagle Summit",]
spec_all = spec_all[,400:850]

names(spec_all) = meta(spec_all)$Species_ID
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 850, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 50)

perf.plot_species = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
                         legend.position = "horizontal")
saveRDS(perf.plot_species, "Figures/perf plots/perf.plot_species.rds")
###ncomp = 16
plotIndiv(plsda.fit, title = "", comp = c(1,2,5), legend = TRUE, 
          style = "3d", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 30)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,19] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)


#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = .75)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)
################################################################################
#Ala and Oct - no hyb
################################################################################

#data
spec_all1 = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_all = spec_all1[!meta(spec_all1)$Species_ID == "DX",]
names(spec_all) = meta(spec_all)$Species_ID
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)
perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 50)

perf.plot_species.nohyb = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
                         legend.position = "horizontal")
par(mfrow = c(1,1))
perf.plot_species.nohyb
saveRDS(perf.plot_species.nohyb, "Figures/perf plots/perf.plot_species.nohyb.rds")
###ncomp = 16
plotIndiv(plsda.fit, title = "", comp = c(1,2,5), legend = TRUE, 
          style = "3d", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 21)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,21] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)


#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)

plotLoadings(plsda.fit, contrib = "max", method = "mean", comp = 1,
             plot = TRUE, show.ties = TRUE, col.ties="white", ndisplay = NULL, size.name = 0.7,
             size.legend = 0.8, name.var=NULL, name.var.complete=FALSE, title = NULL,
             size.title = rel(1.8), size.subtitle = rel(1.4),
             legend = TRUE, legend.color = NULL, legend.title = 'Outcome',
             layout = NULL, border = NA, xlim = c(-.15,.15))


plsda.load = plotLoadings(plsda.fit, contrib = "max", comp = 1, plot = FALSE)
plsda.load = plsda.load[ sort(rownames(plsda.load)), ]
barplot(plsda.load$importance, main = "Loadings for Species", col = plsda.load$color, border = FALSE, 
        lwd = 0.01, space = 0, ylab = "Importance", xlab = "Wavelength (nm)", 
        ylim = c(-.10, .15))
rect(1,-.3,31, .15)
rect(41,-.3,91, .15)
rect(116, -.3, 141, .15)
rect(161,-.3, 201, .15)

axis(2)
axis(1,  at = c(1,11, 61, 111, 161,200), labels = c(400,500, 1000, 1500, 2000,2400))
axis(3, at=c(16,66,128.5,181), tick = FALSE, labels = c("VIS", "NIR", "SWIR 1", "SWIR 2"),
     line = -1)
legend(180, .14, title="Species",
       c("DO", "DA"), fill= c("#F68B33", "#388ECC"), horiz=F, cex=0.8)




x = plsda.fit$loadings$X
Comp1 = x[,1]
Comp2 = x[,2]
Comp3 = x[,3]

par(mar = c(5,5,5,1), mgp = c(3, 1, 0))
plot(Comp1, main = "PLSDA Loadings for Species - Dry, No Hybrids", xaxt = 'n', type = "l", lwd = 3,
     col = "red", xlab = "Wavelength", ylab = "Loading Value", ylim = c(-0.25, 0.20),
     panel.first = c(rect(1,-.3,31, .25, col = "gray95", border = NA),
                     rect(41,-.3,91, .25, col = "gray95", border = NA),
                     rect(116, -.3, 141, .25, col = "gray95", border = NA),
                     rect(161,-.3, 201, .25, col = "gray95", border = NA)))
axis(1, at=c(11, 61, 111, 161), labels=rownames(x)[c(11, 61, 111, 161)])
axis(3, at=c(16,66,128.5,181), tick = FALSE, labels = c("VIS", "NIR", "SWIR 1", "SWIR 2"),
     line = -1)
lines(Comp2, col="turquoise2",lty=1, lwd = 3)
lines(Comp3, col="purple", lty = 1, lwd = 3)
abline(h = 0, lty = 2, lwd = 2)
legend("bottomright", inset=.02, title="Components",
       c("1","2","3"), fill= c("red", "turquoise2", "purple"), horiz=TRUE, cex=0.8)

################################################################################
# Fit PLS_DA model: predict location
################################################################################
#data
spec_all = readRDS("Clean-up/Vector_normalized/all_vn.rds")
names(spec_all) = meta(spec_all)$Location
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 50)

perf.plot_loc = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
saveRDS(perf.plot_loc, "Figures/perf plots/perf.plot_loc.rds")

###ncomp = 26
plotIndiv(plsda.fit, title = "", comp = c(1,16), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 24)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,24] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)


#plot
par(mar = c(2, 5, 3, 1), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)

#Loadings Plots
x = plsda.fit$loadings$X
Comp1 = x[,1]
Comp2 = x[,2]
Comp3 = x[,3]

par(mar = c(5,5,5,1), mgp = c(3, 1, 0))
plot(Comp1, main = "PLSDA Loadings for Location - Dry", xaxt = 'n', type = "l", lwd = 3,
     col = "red", xlab = "Wavelength", ylab = "Loading Value", ylim = c(-0.25, 0.20),
     panel.first = c(rect(1,-.3,31, .25, col = "gray95", border = NA),
                     rect(41,-.3,91, .25, col = "gray95", border = NA),
                     rect(116, -.3, 141, .25, col = "gray95", border = NA),
                     rect(161,-.3, 201, .25, col = "gray95", border = NA)))
axis(1, at=c(11, 61, 111, 161), labels=rownames(x)[c(11, 61, 111, 161)])
axis(3, at=c(16,66,128.5,181), tick = FALSE, labels = c("VIS", "NIR", "SWIR 1", "SWIR 2"),
     line = -1)
lines(Comp2, col="turquoise2",lty=1, lwd = 3)
lines(Comp3, col="purple", lty = 1, lwd = 3)
abline(h = 0, lty = 2, lwd = 2)
legend("bottomright", inset=.02, title="Components",
       c("1","2","3"), fill= c("red", "turquoise2", "purple"), horiz=TRUE, cex=0.8)


abline(v = 106, lty = 1, lwd = 3, col = "blue")
abline(v = 155, lty = 1, lwd = 3, col = "blue")
abline(v = 59, lty = 1, lwd = 1, col = "blue")
abline(v = 85, lty = 1, lwd = 1, col = "blue")
################################################################################
# Fit PLS_DA model all (species plus location)
################################################################################
#data
spec_all = readRDS("lean-up/Vector_normalized/all_vn.rds")
names(spec_all) = meta(spec_all)$Species
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)

#add location as variable
spec_loc_mat_s = cbind(as.factor(spec_all.df$Location), spec_mat_s)
spec_mat   = spec_loc_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 20)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 16
plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)


#Run PLSDA
set.seed(26) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 16)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,16] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)


#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)

################################################################################
# Fit PLS_DA model - predicting location of octopetala
################################################################################
#data
spec_all = readRDS("Clean-up/Vector_normalized/all_vn.rds")
oct = spec_all[meta(spec_all)$Species == "octopetala",]
names(oct) = meta(oct)$Location
oct.df = as.data.frame(oct)

#resample by every 10 nm
spec_small = resample(oct, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)
###ncomp = 24

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 24)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,24] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)


#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = .55)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = .55)


################################################################################
# Fit PLS_DA model 3 sites, 3 species
################################################################################

#data
spec = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec1 = spec[meta(spec)$Location == "Wickersham Dome B",]
spec2 = spec[meta(spec)$Location == "Wickersham Dome A",]
spec3 = spec[meta(spec)$Location == "Eagle Summit",]
big3 = Reduce(combine, list(spec1, spec2, spec3))
names(big3) = meta(big3)$Species_ID
big3.m = as.matrix(big3)
big3.df = as.data.frame(big3)

#Resample by 10 nm
spec_small = resample(big3, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

###Ncomp = 25
plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(27) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 30)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,24] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)

#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)


################################################################################
# Fit PLS_DA model 3 sites, 3 species wet
################################################################################

#data
big3 = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
names(big3) = meta(big3)$Species_ID
big3.m = as.matrix(big3)
big3.df = as.data.frame(big3)

#Resample by 10 nm
spec_small = resample(big3, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

###Ncomp = 25
plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(27) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 30)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,20] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)

#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)
################################################################################
# Fit PLS_DA model 3 sites, 3 species plus location
################################################################################

#data
big3 = readRDS("Clean-up/Vector_normalized/vn_big3.rds")
names(big3) = meta(big3)$Species
big3.m = as.matrix(big3)
big3.df = as.data.frame(big3)

#resample by every 10 nm
spec_small = resample(big3, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)

#add location as variable
spec_loc_mat_s = cbind(as.factor(big3.df$Location), spec_mat_s)
spec_mat = spec_loc_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#Determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)
###Ncomp = 25
#Run PLSDA
set.seed(28) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 25)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,25] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)

#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)

################################################################################
# Fit PLS_DA model 3 sites, 3 species - predict location
################################################################################

#data
big3 = readRDS("Clean-up/Vector_normalized/vn_big3.rds")
names(big3) = meta(big3)$Location
big3.m = as.matrix(big3)
big3.df = as.data.frame(big3)

#Resample by 10 nm
spec_small = resample(big3, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 13)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

###Ncomp = 13
plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(27) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 13)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,25] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)

#plot
par(mar = c(2, 5, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)


################################################################################
# Fit PLS_DA model 3 sites, no hybrids
################################################################################
#Data
no_hybrids = readRDS("Clean-up/Vector_normalized/vn_big3.no_hybrids.rds")
names(no_hybrids) = meta(no_hybrids)$Species_ID
no_hybrids.m = as.matrix(no_hybrids)
no_hybrids.df = as.data.frame(no_hybrids)

#resample by 10 nm
spec_small = resample(no_hybrids, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 50)

saveRDS(perf.plsda, "perf.plsda.rds")

perf.plot_nohyb = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
saveRDS(perf.plot_nohyb, "Figures/perf plots/perf.plot_nohyb.rds")

plotIndiv(plsda.fit, title = "", comp = c(1,2,3), legend = TRUE, 
          style = "3d", ind.names = F, ellipse = TRUE)
###Ncomp = 20
#run plsda
set.seed(2543) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 20)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,20] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)

#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)

################################################################################
# Fit PLS_DA model 3 sites, no hybrids plus location
################################################################################
#data
no_hybrids = readRDS("Clean-up/Vector_normalized/vn_big3.no_hybrids.rds")
names(no_hybrids) = meta(no_hybrids)$Species
no_hybrids.m = as.matrix(no_hybrids)
no_hybrids.df = as.data.frame(no_hybrids)

#Resample by every 10 nm
spec_small = resample(no_hybrids, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_loc_mat_s = cbind(as.factor(no_hybrids.df$Location), spec_mat_s)
spec_mat = spec_loc_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "vertical")
###ncomp = 16
plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

set.seed(29) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 16)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,16] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)

#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)


################################################################################
# Fit PLS_DA model 3 sites, no hybrids - predict location
################################################################################
#Data
no_hybrids = readRDS("Clean-up/Vector_normalized/vn_big3.no_hybrids.rds")
names(no_hybrids) = meta(no_hybrids)$Location
no_hybrids.m = as.matrix(no_hybrids)
no_hybrids.df = as.data.frame(no_hybrids)

#resample by 10 nm
spec_small = resample(no_hybrids, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 30)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

saveRDS(perf.plsda, "perf.plsda.rds")

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,3), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)
###Ncomp = 11
#run plsda
set.seed(2543) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 11)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,11] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)

#plot
par(mar = c(2, 4, 3, 4), oma = c(2, 4, 3, 2))
color2D.matplot(cm1, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)


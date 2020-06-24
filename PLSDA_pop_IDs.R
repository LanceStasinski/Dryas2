#PLSDA for predicting pop ID
#RDS files for cleaned data can be downloaded from:
#FM_Polar_Studies > Dryas_Spectral_Analyses > Clean-up 
#or generated via 01_clean_up.R
################################################################################
#Set up
################################################################################

library(mixOmics)
library(spectrolab)
library(plotrix)

################################################################################
#Set working directory to folder containing downloaded rds files
################################################################################

setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
# Fit PLS_DA model all dry
################################################################################
#data
spec_all1 = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_all = spec_all1[!meta(spec_all1)$GenePop_ID == "NaN",]
names(spec_all) = meta(spec_all)$GenePop_ID
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

perf.plot_pops = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

saveRDS(perf.plot_pops, "Figures/perf plots/perf.plot_pops.rds")
###ncomp = 23
plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 26)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,26] 
# calculate the error rate of the model
confusion.mat = get.confusion_matrix(truth = resp[test], predicted = prediction)
cm1 = as.data.frame(confusion.mat)
get.BER(confusion.mat)


#plot
par(mar = c(2, 4, 3, 1), oma = c(2, 4, 3, 2))
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
             layout = NULL, border = NA, xlim = NULL)


x = plsda.fit$loadings$X
Comp1 = x[,1]
Comp2 = x[,2]
Comp3 = x[,3]

par(mar = c(5,5,5,1), mgp = c(3, 1, 0))
plot(Comp1, main = "PLSDA Loadings for Population ID - Dry", xaxt = 'n', type = "l", lwd = 3,
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
# Fit PLS_DA model all wet
################################################################################
#data
spec_all1 = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec_all = spec_all1[!meta(spec_all1)$GenePop_ID == "NaN",]
names(spec_all) = meta(spec_all)$GenePop_ID
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

perf.plot_pops_wet = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
saveRDS(perf.plot_pops_wet, "Figures/perf plots/perf.plot_pops_wet.rds")
###ncomp = 19
plotIndiv(plsda.fit, title = "", comp = c(4,2,3), legend = TRUE, 
          style = "3d", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(25) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 18)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,18] 
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

#Loadings plot

x = plsda.fit$loadings$X
Comp1 = x[,1]
Comp2 = x[,2]
Comp3 = x[,3]

par(mar = c(5,5,5,1), mgp = c(3, 1, 0))
plot(Comp1, main = "PLSDA Loadings for Population ID - Wet", xaxt = 'n', type = "l", lwd = 3,
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
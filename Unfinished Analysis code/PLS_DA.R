#PCA
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
# Fit PLS_DA model all
################################################################################
#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
names(spec_all) = meta(spec_all)$Species
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
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 18
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

################################################################################
# Fit PLS_DA model: predict location
################################################################################
#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
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
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 18
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
# Fit PLS_DA model all (species plus location)
################################################################################
#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
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
###ncomp = 12
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
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 12)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,12] 
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
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
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
###ncomp = 25

#Run PLSDA
set.seed(25) 
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
     labels = names(cm1), tick = FALSE, cex.axis = .55)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = .55)


################################################################################
# Fit PLS_DA model 3 sites, 3 species
################################################################################

#data
big3 = readRDS("Clean-up/Clean_spectra/clean_big3.rds")
names(big3) = meta(big3)$Species
big3.m = as.matrix(big3)
big3.df = as.data.frame(big3)

#Resample by 10 nm
spec_small = resample(big3, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

#determine number of components to use
plsda.fit = plsda(spec_mat, resp, ncomp = 25)

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
# Fit PLS_DA model 3 sites, 3 species plus location
################################################################################

#data
big3 = readRDS("Clean-up/Clean_spectra/clean_big3.rds")
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
# Fit PLS_DA model 3 sites, no hybrids
################################################################################
#Data
no_hybrids = readRDS("Clean-up/Clean_spectra/big3.no_hybrids.rds")
names(no_hybrids) = meta(no_hybrids)$Species
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
###Ncomp = 15
#run plsda
set.seed(2543) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 

# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 15)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,15] 
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
no_hybrids = readRDS("Clean-up/Clean_spectra/big3.no_hybrids.rds")
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
plsda.fit = plsda(spec_mat, resp, ncomp = 20)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "vertical")
###ncomp = 18
plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

set.seed(29) 
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



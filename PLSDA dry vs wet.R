#PLSDA Wet vs Dry
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
# Fit PLS_DA model wet ES
################################################################################
#data
spec_wet = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec_es = spec_wet[meta(spec_wet)$Location == "Eagle Summit",]
spec_all = spec_es
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
plsda.fit = plsda(spec_mat, resp, ncomp = 40)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 25
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
# Fit PLS_DA model wet ES nohyb
################################################################################
#data
spec_wet = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec_es = spec_wet[meta(spec_wet)$Location == "Eagle Summit",]
spec_es_nohyb = spec_es[!meta(spec_es)$Species == "hybrid",]
spec_all = spec_es_nohyb
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
plsda.fit = plsda(spec_mat, resp, ncomp = 40)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 22
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
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 22)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,22] 
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
#PLSDA ES dry 
################################################################################

#data
spec_dry = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_es_d = spec_dry[meta(spec_dry)$Location == "Eagle Summit",]
spec_all = spec_es_d
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
plsda.fit = plsda(spec_mat, resp, ncomp = 40)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 24
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
     labels = names(cm1), tick = FALSE, cex.axis = 1)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)
################################################################################
# Fit PLS_DA model dry ES no hyb
################################################################################
#data
spec_dry = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_es_d = spec_dry[meta(spec_dry)$Location == "Eagle Summit",]
spec_es_d_nohyb = spec_es_d[!meta(spec_es_d)$Species == "hybrid",]
spec_all = spec_es_d_nohyb
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
plsda.fit = plsda(spec_mat, resp, ncomp = 40)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 15
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
#PLSDA wdb wet
################################################################################
#data
spec_wet = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec_es = spec_wet[meta(spec_wet)$Location == "Wickersham Dome B",]
spec_all = spec_es
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
plsda.fit = plsda(spec_mat, resp, ncomp = 40)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 31
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
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 31)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,31] 
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
# Fit PLS_DA model wet wdb nohyb
################################################################################
#data
spec_wet = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec_wdb = spec_wet[meta(spec_wet)$Location == "Wickersham Dome B",]
spec_wdb_nohyb = spec_wdb[!meta(spec_wdb)$Species == "hybrid",]
spec_all = spec_wdb_nohyb
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
plsda.fit = plsda(spec_mat, resp, ncomp = 40)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 14
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
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 14)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,14] 
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
#PLSDA wdb dry 
################################################################################

#data
spec_dry = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_wdb_d = spec_dry[meta(spec_dry)$Location == "Wickersham Dome B",]
spec_all = spec_wdb_d
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
plsda.fit = plsda(spec_mat, resp, ncomp = 40)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 25
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
# Fit PLS_DA model dry wdb no hyb
################################################################################
#data
spec_dry = readRDS("Clean-up/Vector_normalized/all_vn.rds")
spec_wdb_d = spec_dry[meta(spec_dry)$Location == "Wickersham Dome B",]
spec_wdb_d_nohyb = spec_wdb_d[!meta(spec_wdb_d)$Species == "hybrid",]
spec_all = spec_wdb_d_nohyb
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
plsda.fit = plsda(spec_mat, resp, ncomp = 40)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")
###ncomp = 21
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
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 21)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the nth component
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
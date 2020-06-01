################################################################################
#Set up
################################################################################

library(mixOmics)
library(spectrolab)
library(plotrix)

#all
spec_all = readRDS("clean_all.rds")
names(spec_all) = meta(spec_all)$Species
spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#octotpetala
spec_all = readRDS("clean_all.rds")
oct = spec_all[meta(spec_all)$Species == "octopetala",]
names(oct) = meta(oct)$Location
oct.df = as.data.frame(oct)

#3 sites, all species
big3 = readRDS("clean_big3.rds")
names(big3) = meta(big3)$Species
big3.m = as.matrix(big3)
big3.df = as.data.frame(big3)

#3 sites, no hybrids
no_hybrids = readRDS("big3.no_hybrids.rds")
names(no_hybrids) = meta(no_hybrids)$Species
no_hybrids.m = as.matrix(no_hybrids)
no_hybrids.df = as.data.frame(no_hybrids)

################################################################################
# Fit PLS_DA model all
################################################################################
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))


plsda.fit = plsda(spec_mat, resp, ncomp = 15)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

set.seed(25) 
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

################################################################################
# Fit PLS_DA model all (species plus location)
################################################################################

spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)

spec_loc_mat_s = cbind(as.factor(spec_all.df$Location), spec_mat_s)

spec_mat   = spec_loc_mat_s


resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

plsda.fit = plsda(spec_mat, resp, ncomp = 20)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

set.seed(26) 
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

################################################################################
# Fit PLS_DA model - predicting location
################################################################################
spec_small = resample(oct, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s
resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))


plsda.fit = plsda(spec_mat, resp, ncomp = 20)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

set.seed(25) 
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
     labels = names(cm1), tick = FALSE, cex.axis = .55)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = .55)


################################################################################
# Fit PLS_DA model 3 sites, 3 species
################################################################################
spec_small = resample(big3, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s

resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

plsda.fit = plsda(spec_mat, resp, ncomp = 20)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

set.seed(27) 
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
# Fit PLS_DA model 3 sites, 3 species plus location
################################################################################
spec_small = resample(big3, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_loc_mat_s = cbind(as.factor(big3.df$Location), spec_mat_s)

spec_mat = spec_loc_mat_s

resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

plsda.fit = plsda(spec_mat, resp, ncomp = 20)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

set.seed(28) 
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
# Fit PLS_DA model 3 sites, equal
################################################################################

spec_small = resample(no_hybrids, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s

resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

plsda.fit = plsda(spec_mat, resp, ncomp = 15)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

saveRDS(perf.plsda, "perf.plsda.rds")

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "horizontal")

plotIndiv(plsda.fit, title = "", comp = c(1,3), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)


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
# Fit PLS_DA model 3 sites, equal plus location
################################################################################
spec_small = resample(no_hybrids, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_loc_mat_s = cbind(as.factor(no_hybrids.df$Location), spec_mat_s)

spec_mat = spec_loc_mat_s

resp = rownames(spec_mat)
rownames(spec_mat) = seq(nrow(spec_mat))

plsda.fit = plsda(spec_mat, resp, ncomp = 20)

perf.plsda = perf(plsda.fit, validation = "Mfold", folds = 5,
                  progressBar = TRUE, auc = TRUE, nrepeat = 10)

plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
     legend.position = "vertical")

plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

set.seed(29) 
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



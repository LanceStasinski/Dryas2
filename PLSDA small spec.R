

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
# PLSDA simulate airborne sensor
################################################################################

#data
spec_all = readRDS("Clean-up/Vector_normalized/vn_all_w.rds")
spec_all = spec_all[meta(spec_all)$Location == "Wickersham Dome B",]
spec_all = spec_all[!meta(spec_all)$Species_ID == "DX",]
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
                extremes = c("white", "deepskyblue3"),
               )
axis(3, at = seq_len(ncol(cm1)) - 0.5,
     labels = names(cm1), tick = FALSE, cex.axis = .75, line = -.5)
axis(2, at = seq_len(nrow(cm1)) -0.5,
     labels = rev(rownames(cm1)), tick = FALSE, las = 1, cex.axis = 1)
mtext("Species prediction for WDB using 400:850 nm (20 comps, 0 BER), No Hybrids",
      side = 3, line = 2)

#Predict species + location
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
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
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

perf.plot_pops = plot(perf.plsda, col = color.mixo(1:3), sd = TRUE, 
                      legend.position = "horizontal")


###ncomp = 23
plotIndiv(plsda.fit, title = "", comp = c(1,2), legend = TRUE, 
          style = "graphics", ind.names = F, ellipse = TRUE)

#Run PLSDA
set.seed(10) 
samp <- sample(1:3, nrow(spec_mat), replace = TRUE) 
# 1/3 of the data will compose the test set
test <- which(samp == 1) 
# rest will compose the training set
train <- setdiff(1:nrow(spec_mat), test)

## For PLS-DA, train the model
plsda.train <- plsda(spec_mat[train, ], resp[train], ncomp = 27)
# then predict
test.predict <- predict(plsda.train, spec_mat[test, ], dist = "max.dist")
# store prediction for the 4th component
prediction <- test.predict$class$max.dist[,27] 
# calculate the error rate of the model
mat10 = get.confusion_matrix(truth = resp[test], predicted = prediction)

mat.total = (mat1 + mat2 + mat3 + mat4 + mat5 + mat6 + mat7 + mat8 + mat9 + mat10)/10
saveRDS(mat.total, "sploc_nohyb_mat.rds")
pops.cm = as.data.frame(mat.total)
get.BER(mat.total)

names(pops.cm)


#plot
par(mar = c(2, 5, 5, 1), oma = c(2, 4, 3, 2))
color2D.matplot(pops.cm, 
                show.values = TRUE,
                axes = FALSE,
                xlab = "",
                ylab = "",
                vcex = 2,
                vcol = "black",
                extremes = c("white", "deepskyblue3"))
axis(3, at = seq_len(ncol(pops.cm)) - 0.5,
     labels = c("DA_es", "DA_tm", "DA_wdb", "DO_bg", "DO_es", "DO_mdb", "DO_tm", "DO_wda", "DO_wdb"),
     tick = FALSE, cex.axis = 1, line = -.5)
axis(2, at = seq_len(nrow(pops.cm)) -0.5,
     labels = rev(rownames(pops.cm)), tick = FALSE, las = 1, cex.axis = 1, line = -.5)
mtext("Ncomps:27    BER:0.02687525    Iterations:10", side = 1,
      line = 1)
mtext("Actual Species and Location", side = 2, line = 5, cex = 1.5)
mtext("Predicted Species and Location", side = 3, line = 3, cex =1.5)
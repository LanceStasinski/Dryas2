#is50 PLSDA
################################################################################
#Set up
################################################################################

library(spectrolab)
library(caret)
library(mlbench)
library(corrplot)
library(matrixStats)
library(naniar)


setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
ispec = readRDS('Clean-up/Clean_spectra/spec_iS50.rds')
spec_all = ispec

spec_all = spec_all[!meta(spec_all)$GenePop_ID == "NaN"]
spec_all = spec_all[!meta(spec_all)$Species_ID == "DX"]


spec_all.m = as.matrix(spec_all)
spec_mat = spec_all.m
spec_all.df = as.data.frame(spec_all)


#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$Species_ID)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$Species_ID"] <- "Species_ID"


for(i in 1:10){
  
  set.seed(i)
  
  inTrain <- caret::createDataPartition(
    y = spec_df$Species_ID,
    p = .8,
    list = FALSE
  )
  
  training <- spec_df[inTrain,]
  testing <- spec_df[-inTrain,]
  
  #tune model
  ctrl <- trainControl(
    method = "repeatedcv", 
    number = 10,
    repeats = 3)
  
  
  plsFit <- train(
    Species_ID ~ .,
    data = training,
    method = "pls",
    preProc = c("center", "scale"),
    trControl = ctrl,
    tuneLength = 24)
  
  assign(paste0('plsFit', i), plsFit)
  
  loadings = plsFit$finalModel$loadings
  loadings.m = as.matrix(loadings)
  class(loadings.m) <- 'matrix'
  assign(paste0('lm',i), loadings.m)
  
  comp1 = loadings.m[,1]
  assign(paste0('comp1_',i), comp1)
  comp2 = loadings.m[,2]
  assign(paste0('comp2_',i), comp2)
  comp3 = loadings.m[,3]
  assign(paste0('comp3_',i), comp3)
  
  
  
  #test model
  plsClasses <- predict(plsFit, newdata = testing)
  
  
  #Confusion matrices
  cm = confusionMatrix(data = plsClasses, testing$Species_ID)
  acc = cm$overall[1]
  assign(paste0('acc',i), acc)
  
  
  cm.m = as.matrix(cm)
  
  assign(paste0("cm", i), cm.m)
}

acc = c(acc1, acc2, acc3, acc4, acc5, acc6, acc7, acc8, acc9, acc10)
mean.acc = mean(acc)
sd.acc = sd(acc)

mean.acc
sd.acc

#Kappa
k1 = as.matrix(plsFit1$results$Kappa)
k2 = as.matrix(plsFit2$results$Kappa)
k3 = as.matrix(plsFit3$results$Kappa)
k4 = as.matrix(plsFit4$results$Kappa)
k5 = as.matrix(plsFit5$results$Kappa)
k6 = as.matrix(plsFit6$results$Kappa)
k7 = as.matrix(plsFit7$results$Kappa)
k8 = as.matrix(plsFit8$results$Kappa)
k9 = as.matrix(plsFit9$results$Kappa)
k10 = as.matrix(plsFit10$results$Kappa)
k.total = Reduce(cbind, list(k1,k2,k3,k4,k5,k6,k7,k8,k9,k10))

kavg = as.matrix(rowMeans(k.total))
ksd = as.matrix(rowSds(k.total))
klower = kavg - ksd
khigher = kavg + ksd

plot(kavg, type = 'p', pch = 16, cex = .75, ylab = 'Kappa', xlab = 'Component', 
     xlim = c(0,100), main = 'Kappa for Species_ID')
lines(klower, lty = 2, col = 'red')
lines(khigher, lty = 2, col = 'red')
abline(v = 19, col = 'blue')
legend('bottomright', legend = c('Mean', 'Standard deviation', 'Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 2, 1), col = c('black', 'red', 'blue'))
#accuracy
a1 = as.matrix(plsFit1$results$Accuracy)
a2 = as.matrix(plsFit2$results$Accuracy)
a3 = as.matrix(plsFit3$results$Accuracy)
a4 = as.matrix(plsFit4$results$Accuracy)
a5 = as.matrix(plsFit5$results$Accuracy)
a6 = as.matrix(plsFit6$results$Accuracy)
a7 = as.matrix(plsFit7$results$Accuracy)
a8 = as.matrix(plsFit8$results$Accuracy)
a9 = as.matrix(plsFit9$results$Accuracy)
a10 = as.matrix(plsFit10$results$Accuracy)
a.total = Reduce(cbind, list(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10))

a.avg = as.matrix(rowMeans(a.total))
a.sd = as.matrix(rowSds(a.total))
alower = a.avg - a.sd
ahigher = a.avg + a.sd

plot(a.avg, type = 'p', pch = 16, cex = .75, ylim = c(.6, 1), ylab = 'Accuracy', xlab = 'Component', 
     xlim = c(0,40), main = 'Accuracy for Species_ID')
lines(alower, lty = 2, col = 'red')
lines(ahigher, lty = 2, col = 'red')
abline(v = 19, col = 'blue')
legend('bottomright', legend = c('Mean', 'Standard deviation', 'Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 2, 1), col = c('black', 'red', 'blue'))



#Confusion matrix

cm.total = (cm1 + cm2 + cm3 + cm4 + cm5 + cm6 + cm7 + cm8+ cm9 + cm10)/10
cm.total = as.matrix(cm)
cm.total = t(cm.total)
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('DA', 'DO_et', 'DO_wdb')


par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
corrplot(cm.total, is.corr = F, method = 'color', addCoef.col = 'darkorange2',
         tl.srt = 0, tl.offset = 1.5, number.digits = 2, tl.cex = .75,
         tl.col = 'black', cl.pos = 'n', na.label = 'square', 
         na.label.col = 'white', addgrid.col = 'grey')
mtext("Reference", side = 2, line = -5, cex = 1.5)
mtext("Prediction", side = 3, cex = 1.5, at = 2, line = 1)

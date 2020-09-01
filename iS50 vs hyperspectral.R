library(spectrolab)
library(caret)
library(mlbench)
library(corrplot)
library(matrixStats)
library(naniar)


setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")
ispec = readRDS('Clean-up/Clean_spectra/spec_iS50.rds')
dry_spec = readRDS('Clean-up/Clean_spectra/clean_all.rds')

names = unique(meta(ispec)$Name)

s1 = dry_spec[meta(dry_spec)$Name == 'DryalaESA20',]

s3 = dry_spec[meta(dry_spec)$Name == 'DryalaESA22',]
s4 = dry_spec[meta(dry_spec)$Name == 'DryalaESA23',]

s6 = dry_spec[meta(dry_spec)$Name == 'DryalaESA25',]
s7 = dry_spec[meta(dry_spec)$Name == 'DryalaESA26',]
s8 = dry_spec[meta(dry_spec)$Name == 'DryalaESA27',]

s10 = dry_spec[meta(dry_spec)$Name == 'DryalaESA29',]
s11 = dry_spec[meta(dry_spec)$Name == 'DryoctESA50',]
s12 = dry_spec[meta(dry_spec)$Name == 'DryoctESA51',]
s13 = dry_spec[meta(dry_spec)$Name == 'DryoctESA52',]

s15 = dry_spec[meta(dry_spec)$Name == 'DryoctESA54',]

s17 = dry_spec[meta(dry_spec)$Name == 'DryoctESA56',]

s19 = dry_spec[meta(dry_spec)$Name == 'DryoctESA58',]
s20 = dry_spec[meta(dry_spec)$Name == 'DryoctESA59',]
s21 = dry_spec[meta(dry_spec)$Name == 'DryalaTM11',]

s23 = dry_spec[meta(dry_spec)$Name == 'DryalaTM13',]
s24 = dry_spec[meta(dry_spec)$Name == 'DryalaTM14',]
s25 = dry_spec[meta(dry_spec)$Name == 'DryalaTM15',]
s26 = dry_spec[meta(dry_spec)$Name == 'DryalaTM16',]
s27 = dry_spec[meta(dry_spec)$Name == 'DryalaTM17',]
s28 = dry_spec[meta(dry_spec)$Name == 'DryalaTM18',]
s29 = dry_spec[meta(dry_spec)$Name == 'DryalaTM19',]
s30 = dry_spec[meta(dry_spec)$Name == 'DryalaTM20',]
s31 = dry_spec[meta(dry_spec)$Name == 'DryoctTM21',]
s32 = dry_spec[meta(dry_spec)$Name == 'DryoctTM22',]
s33 = dry_spec[meta(dry_spec)$Name == 'DryoctTM23',]
s34 = dry_spec[meta(dry_spec)$Name == 'DryoctTM24',]
s35 = dry_spec[meta(dry_spec)$Name == 'DryoctTM25',]
s36 = dry_spec[meta(dry_spec)$Name == 'DryoctTM26',]
s37 = dry_spec[meta(dry_spec)$Name == 'DryoctTM27',]
s38 = dry_spec[meta(dry_spec)$Name == 'DryoctTM28',]
s39 = dry_spec[meta(dry_spec)$Name == 'DryoctTM29',]
s40 = dry_spec[meta(dry_spec)$Name == 'DryoctTM30',]
s41 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB1',]
s42 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB2',]
s43 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB3',]
s44 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB4',]
s45 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB5',]
s46 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB6',]
s47 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB7',]
s48 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB8',]
s49 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB9',]
s50 = dry_spec[meta(dry_spec)$Name == 'DryalaWDB10',]
s51 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB21',]
s52 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB22',]
s53 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB23',]
s54 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB24',]
s55 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB25',]
s56 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB26',]
s57 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB28',]
s58 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB29',]
s59 = dry_spec[meta(dry_spec)$Name == 'DryoctWDB30',]

dry3 = Reduce(combine, list(s1,s3,s4,s6,s7,s8,s10,s12,s13,s15,s17,s19,s20,s21,s23,
                            s24,s25,s26,s27,s28,s29,s30,s31,s32,s33,s34,s35,s36,
                            s37,s38,s39,s40,s41,s42,s43,s44,s45,s46,s47,s48,s49,
                            s50,s51,s52,s53,s54,s55,s56,s57,s58,s59))

ispec = ispec[!meta(ispec)$Name == "DryalaESA21",]
ispec = ispec[!meta(ispec)$Name == "DryalaESA24",]
ispec = ispec[!meta(ispec)$Name == "DryalaESA28",]
ispec = ispec[!meta(ispec)$Name == "DryoctESA53",]
ispec = ispec[!meta(ispec)$Name == "DryoctESA55",]
ispec = ispec[!meta(ispec)$Name == "DryoctESA57",]
ispec = ispec[!meta(ispec)$Name == "DryalaTM12",]




library(spectrolab)
library(caret)
library(mlbench)
library(corrplot)
library(matrixStats)
library(naniar)


setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
# Fit PLS_DA model all dry
################################################################################

#data
spec_all = dry3

spec_all = ispec

spec_all.m = as.matrix(spec_all)
spec_mat = spec_
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s

#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$Location)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$Location"] <- "Location"


#Partition Data

for(i in 1:10){
  
  set.seed(i)
  
  inTrain <- caret::createDataPartition(
    y = spec_df$Location,
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
    Location ~ .,
    data = spec_df,
    method = "pls",
    preProc = c("center", "scale"),
    trControl = ctrl,
    tuneLength = 60)
  
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
  cm = confusionMatrix(data = plsClasses, testing$Location)
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



#kappa

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
     xlim = c(0,40), main = 'Kappa for Location')
lines(klower, lty = 2, col = 'red')
lines(khigher, lty = 2, col = 'red')
abline(v = 25, col = 'blue')
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

plot(a.avg, type = 'p', pch = 16, cex = .75, ylab = 'Accuracy', xlab = 'Component', 
     xlim = c(0,40), main = 'Accuracy for Location')
lines(alower, lty = 2, col = 'red')
lines(ahigher, lty = 2, col = 'red')
abline(v = 16, col = 'blue')
legend('bottomright', legend = c('Mean', 'Standard deviation', 'Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 2, 1), col = c('black', 'red', 'blue'))


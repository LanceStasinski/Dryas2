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
library(rlist)
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")


################################################################################
# set iS50 and PSR+ spectra roughly equal
################################################################################

ispec = readRDS('Clean-up/Clean_spectra/spec_iS50.rds')
hspec = readRDS("Clean-up/Clean_spectra/clean_all.rds")

#Constrain data sets to only include the same samples
ispec.names = meta(ispec)$Name
hspec.names = meta(hspec)$Name
remove.in.ispec = setdiff(ispec.names, hspec.names)
remove.in.hspec = setdiff(hspec.names, ispec.names)
ispec2 = ispec[!meta(ispec)$Name %in% remove.in.ispec,]
hspec2 = hspec[!meta(hspec)$Name %in% remove.in.hspec,]

#Constrain spectra to same range
ispec2 = ispec2[,bands(ispec2, 833, 2400)]
hspec2 = hspec2[,bands(hspec2, 833, 2400)]

################################################################################
# Fit PLS_DA model all dry
################################################################################
spec_all = hspec2

spec_all = spec_all[!meta(spec_all)$GenePop_ID == "NaN"]
spec_all = spec_all[!meta(spec_all)$Species_ID == "DX"]

spec_mat = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$GenePop_ID)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$GenePop_ID"] <- "GenePop_ID"

#Set number of components to be used
ncomp = 16

#Partition Data
accuracy <- c()
kappa <- c()
k.fit <- matrix(nrow = ncomp)
cm.list <- list()

for(i in 1:100){
  
  inTrain <- caret::createDataPartition(
    y = spec_df$GenePop_ID,
    p = .7,
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
    GenePop_ID ~ .,
    data = training,
    maxit = 10000,
    method = "pls",
    trControl = ctrl,
    tuneLength = ncomp)
  
  
  
  #objects for determining n components
  k = assign(paste0('k', i), as.matrix(plsFit$results$Kappa))
  k.fit <- cbind(k.fit, get('k'))
  
  #test model
  plsClasses <- predict(plsFit, newdata = testing)
  
  #objects to assess accuracy 
  cm = confusionMatrix(data = plsClasses, as.factor(testing$GenePop_ID))
  cm.m = assign(paste0("cm", i), as.matrix(cm))
  cm.list <- list.append(cm.list, get('cm.m'))
  
  ac <- assign(paste0('acc',i), cm$overall[1])
  accuracy <- append(accuracy, get('ac'))
  
  kap = assign(paste0("kap",i), cm$overall[2])
  kappa <- append(kappa, get('kap'))
  
}

#accuracy and kappa
mean.acc = mean(accuracy)
sd.acc = sd(accuracy)

mean.kap = mean(kappa)
sd.kap = sd(kappa)

mean.acc
sd.acc

mean.kap
sd.kap

#kappa (for choosing components)

k.total = k.fit[,-1]
kavg = as.matrix(rowMeans(k.total))
ksd = as.matrix(rowSds(k.total))

klower = kavg - ksd
khigher = kavg + ksd
x = 1:30
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, kavg, type = 'p', pch = 16, cex = .75, ylab = 'Kappa', xlab = 'Component', 
     xlim = c(1,60), main = 'Kappa for GenePop_ID')
arrows(x, klower, x, khigher,length=0.05, angle=90, code=3)
abline(v = 16, col = 'blue')
abline(h = max(klower), col = "Red")
legend('bottomright', legend = c('Mean', 'Maximum kappa','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))


#plot confusion matrix
cm.total = Reduce('+', cm.list)/100
cm.total = t(cm.total)
cm.total = cm.total/rowSums(cm.total)

cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('ES', 'TM', 'WDB')
colnames(cm.total) <- c('ES', 'TM', 'WDB')



#plot
pdf(file= "Figures/cm_final/dry/test.pdf", width = 6, height = 6)

dev.new(width = 6, height = 8, unit = 'in')

par(mar = c(1,2,4,1), oma = c(1,1,3,1))
corrplot(cm.total, is.corr = T, method = 'square', addCoef.col = 'darkorange2',
         tl.srt = 0, tl.offset = 1, number.digits = 4, tl.cex = 1.5, 
         cl.cex = 1.5,
         tl.col = 'black', cl.pos = 'n', na.label = 'square', 
         na.label.col = 'white', addgrid.col = 'grey')
mtext("Reference", side = 2, line = -11, cex = 2.5)
mtext("Prediction", side = 3, cex = 2.5, at = 2, line = 5)



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

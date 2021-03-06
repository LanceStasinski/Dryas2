#Caret PLSDA
################################################################################
#Set up
################################################################################

library(spectrolab)
library(caret)
library(dplyr)
library(mlbench)
library(corrplot)
library(matrixStats)
library(naniar)
library(rlist)
setwd("C:/Users/istas/OneDrive/Documents/Dryas Research/Dryas 2.0")

################################################################################
#Data setup !!!!NOTE: use ctrl+f to find a replace the field to be classified!!!
################################################################################

#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all_6scans.rds")

#remove any NaN values and hybrids
spec_all = spec_all[!meta(spec_all)$sp_loc == "NaN",]
spec_all = spec_all[!meta(spec_all)$Species_ID == "DX",]

spec_mat = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#combine relevant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$sp_loc)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$sp_loc"] <- "sp_loc"

################################################################################
#Run PLSDA
################################################################################

#Set number of components to be used
ncomp = 55

#create vectors, lists, and matrices to store metrics and loadings
accuracy <- c()
kappa <- c()
a.fit <- matrix(nrow = ncomp)
cm.list <- list()
da_es.vip = matrix(nrow=2001)
da_tm.vip = matrix(nrow=2001)
da_wdb.vip = matrix(nrow=2001)
do_bg.vip = matrix(nrow=2001)
do_es.vip = matrix(nrow=2001)
do_md.vip = matrix(nrow=2001)
do_tm.vip = matrix(nrow=2001)
do_wda.vip = matrix(nrow=2001)
do_wdb.vip = matrix(nrow=2001)


#start of PLSDA code
for(i in 1:100){
  
  #create data partition: 70% of data for training, 30% for testing
  inTrain <- caret::createDataPartition(
    y = spec_df$sp_loc,
    p = .8,
    list = FALSE
  )
  
  training <- spec_df[inTrain,]
  testing <- spec_df[-inTrain,]
  
  #tune model: 10-fold cross-validation repeated 3 times
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    sampling = 'up',
    repeats = 3)
  
  #Fit model. Note max iterations set to 10000 to allow model convergence
  plsFit <- train(
    sp_loc ~ .,
    data = training,
    maxit = 100000,
    method = "pls",
    trControl = ctrl,
    tuneLength = ncomp)
  
  #variable importance
  vip = varImp(plsFit)
  
  da_es = assign(paste0('da_es', i), vip$importance$DA_es)
  da_es.vip <- cbind(da_es.vip, get('da_es'))
  
  da_tm = assign(paste0('da_tm', i), vip$importance$DA_tm)
  da_tm.vip <- cbind(da_tm.vip, get('da_tm'))
  
  da_wdb = assign(paste0('da_wdb', i), vip$importance$DA_wdb)
  da_wdb.vip <- cbind(da_wdb.vip, get('da_wdb'))
  
  do_bg = assign(paste0('do_bg', i), vip$importance$DO_bg)
  do_bg.vip <- cbind(do_bg.vip, get('do_bg'))
  
  do_es = assign(paste0('do_es', i), vip$importance$DO_es)
  do_es.vip <- cbind(do_es.vip, get('do_es'))
  
  do_md = assign(paste0('do_md', i), vip$importance$DO_mdb)
  do_md.vip <- cbind(do_md.vip, get('do_md'))
  
  do_tm = assign(paste0('do_tm', i), vip$importance$DO_tm)
  do_tm.vip <- cbind(do_tm.vip, get('do_tm'))
  
  do_wda = assign(paste0('do_wda', i), vip$importance$DO_wda)
  do_wda.vip <- cbind(do_wda.vip, get('do_wda'))
  
  do_wdb = assign(paste0('do_wdb', i), vip$importance$DO_wdb)
  do_wdb.vip <- cbind(do_wdb.vip, get('do_wdb'))

  #accuracy objects for determining n components
  a = assign(paste0('a', i), as.matrix(plsFit$results$Accuracy))
  a.fit <- cbind(a.fit, get('a'))
  
  #test model using the testing data partition (30% of data)
  plsClasses <- predict(plsFit, newdata = testing)
  
  #confusion/classification matrix objects to assess accuracy 
  cm = confusionMatrix(data = plsClasses, as.factor(testing$sp_loc))
  cm.m = assign(paste0("cm", i), as.matrix(cm))
  cm.list <- list.append(cm.list, get('cm.m'))
  
  ac <- assign(paste0('acc',i), cm$overall[1])
  accuracy <- append(accuracy, get('ac'))
  
  kap = assign(paste0("kap",i), cm$overall[2])
  kappa <- append(kappa, get('kap'))
}

################################################################################
#Kappa and Accuracy assessed after 100 iterations
################################################################################
mean.acc = mean(accuracy)
sd.acc = sd(accuracy)

mean.kap = mean(kappa)
sd.kap = sd(kappa)

mean.acc
sd.acc

mean.kap
sd.kap

################################################################################
#accuracy values for choosing the optimal number of components to use
################################################################################

a.total = a.fit[,-1]
a.avg = as.matrix(rowMeans(a.total))
a.sd = as.matrix(rowSds(a.total))

a.lower = a.avg - a.sd
a.higher = a.avg + a.sd

#Graph to visually choose optimal number of components
x = 1:60
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, a.avg, type = 'p', pch = 16, cex = .75, ylab = 'Accuracy', 
     xlab = 'Component', xlim = c(1,60), main = 'Accuracy for Species_ID', 
     ylim = c(0,1))
arrows(x, a.lower, x, a.higher,length=0.05, angle=90, code=3)
abline(v = which.max(a.avg), col = 'blue')
abline(h = max(a.avg), col = "Red")
legend('bottomright', legend = c('Mean', 'Maximum accuracy','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))

################################################################################
#Confusion/Classification Matrices
################################################################################
#take average of 100 confusion matrices, reorient matrix, change averages to 
#proportions
cm.avg = Reduce('+', cm.list)/100
cm.avg = t(cm.avg)
cm.total = cm.avg/rowSums(cm.avg)

#standard deviations
f1 <- function(lst){
  n <- length(lst); 	   
  rc <- dim(lst[[1]]); 	   
  ar1 <- array(unlist(lst), c(rc, n)); 	   
  round(apply(ar1, c(1, 2), sd), 2); 	         
}
cm.sd = f1(cm.list)
cm.sd = t(cm.sd)
cm.sd = cm.sd/rowSums(cm.avg)
rownames(cm.sd) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'MD', 
                     'TM', 'WDA', 'WDB')
colnames(cm.sd) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'MD', 
                     'TM', 'WDA', 'WDB')
write.csv(cm.sd, file = 'Figures/cm_final/6_scans/standard deviations/sp_loc_sd_nohyb.csv')

#format matrix for plotting
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'MD', 
                        'TM', 'WDA', 'WDB')
colnames(cm.total) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'MD', 
                        'TM', 'WDA', 'WDB')

#save confusion matrix
write.csv(cm.total, "Figures/cm_final/6_scans/sp_loc_nohyb.csv")

#species + sp_loc special code
cm.total = read.csv("Figures/cm_final/6_scans/sp_loc_nohyb.csv", stringsAsFactors = T)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- cm.total[,1]
cm.total = cm.total[,-1]
cm.total = mapply(cm.total, FUN = as.numeric)
cm.total = matrix(data = cm.total, ncol = 9, nrow = 9)
rownames(cm.total) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'TM', 
                        'MD', 'WDA', 'WDB')
colnames(cm.total) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'TM', 
                        'MD', 'WDA', 'WDB')

#plot confusion matrix

cols = colorRampPalette(c('#f5f5f5', '#fe9929'))

par(mfrow = c(1,1), mar = c(1,2,2,1), oma = c(1,1,3,1))
corrplot::corrplot(cm.total,
                   is.corr = T, 
                   method = 'square', 
                   col = cols(10),
                   addCoef.col = '#542788',
                   tl.srt = 0, 
                   tl.offset = 1, 
                   number.digits = 2, 
                   tl.cex = 1.2, 
                   cl.cex = 1, 
                   number.cex = 1.5,
                   tl.col = 'black', 
                   cl.pos = 'n', 
                   na.label = 'square', 
                   na.label.col = 'white',
                   addgrid.col = 'grey')
mtext("Reference", side = 2, line = -8, cex = 2.5)
mtext("Prediction", side = 3, cex = 2.5, at = 2, line = 3)

################################################################################
#Variable importance
################################################################################
vip_to_spec = function(x){
  t.vip = t(x)
  colnames(t.vip) <- seq(400,2400, by = 1)
  s.vip = as_spectra(t.vip)
}
da_es.vip = da_es.vip[,-1]
da_es.vip.spec = vip_to_spec(da_es.vip)

da_tm.vip = da_tm.vip[,-1]
da_tm.vip.spec = vip_to_spec(da_tm.vip)

da_wdb.vip = da_wdb.vip[,-1]
da_wdb.vip.spec = vip_to_spec(da_wdb.vip)

do_bg.vip = do_bg.vip[,-1]
do_bg.vip.spec = vip_to_spec(do_bg.vip)

do_es.vip = do_es.vip[,-1]
do_es.vip.spec = vip_to_spec(do_es.vip)

do_md.vip = do_md.vip[,-1]
do_md.vip.spec = vip_to_spec(do_md.vip)

do_tm.vip = do_tm.vip[,-1]
do_tm.vip.spec = vip_to_spec(do_tm.vip)

do_wda.vip = do_wda.vip[,-1]
do_wda.vip.spec = vip_to_spec(do_wda.vip)

do_wdb.vip = do_wdb.vip[,-1]
do_wdb.vip.spec = vip_to_spec(do_wdb.vip)


#plot
par(mfrow = c(3,3))
plot(mean(da_es.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = NA, main = 'DAK at ES')
plot_quantile(da_es.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)


plot(mean(da_tm.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = NA, main = 'DAK at TM')
plot_quantile(da_tm.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)



plot(mean(da_wdb.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = NA, main = 'DAK at WDB')
plot_quantile(da_wdb.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)


plot(mean(do_bg.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = NA, main = 'DAJ at BG')
plot_quantile(do_bg.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)




plot(mean(do_es.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = NA, main = 'DAJ at ES')
plot_quantile(do_es.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)



plot(mean(do_md.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA,
     xlab = NA, main = 'DAJ at MD')
plot_quantile(do_md.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)



plot(mean(do_tm.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = 'Wavelength (nm)', main = 'DAJ at TM')
plot_quantile(do_tm.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)


plot(mean(do_wda.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = 'Wavelength (nm)', main = 'DAJ at WDA')
plot_quantile(do_wda.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)



plot(mean(do_wdb.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = 'Wavelength (nm)', main = 'DAJ at WDB')
plot_quantile(do_wdb.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)




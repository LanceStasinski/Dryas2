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
spec_all = normalize(spec_all)
spec_all = spec_all[, 750:2400]
spec_all = spec_all[, -seq(1300, 1460, by = 1)]
spec_all = spec_all[, -seq(1750, 2030, by = 1)]
spec_all = spec_all[, -seq(960, 980, by= 1)]
spec_all = spec_all[, -seq(1170, 1190, by= 1)]
spec_all = spec_all[, -seq(1235, 1255, by= 1)]
spec_all = spec_all[, -seq(2040, 2060, by= 1)]
spec_all = spec_all[, -seq(2135, 2155, by= 1)]
spec_all = spec_all[, -seq(2153, 2173, by= 1)]






#remove any NaN values - mostly pertains to populations
spec_all = spec_all[!meta(spec_all)$Location == "NaN",]

spec_mat = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$Location)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$Location"] <- "Location"


################################################################################
#Run PLSDA
################################################################################

#Set number of components to be used
ncomp = 60

#create vectors, lists, and matrices to store metrics and loadings
accuracy <- c()
kappa <- c()
a.fit <- matrix(nrow = ncomp)
cm.list <- list()
bg.vip = matrix(nrow=1093)
es.vip = matrix(nrow=1093)
md.vip = matrix(nrow=1093)
tm.vip = matrix(nrow=1093)
wda.vip = matrix(nrow=1093)
wdb.vip = matrix(nrow=1093)

#start of PLSDA code
for(i in 1:10){
  
  #create data partition: 70% of data for training, 30% for testing
  inTrain <- caret::createDataPartition(
    y = spec_df$Location,
    p = .8,
    list = FALSE
  )
  
  training <- spec_df[inTrain,]
  testing <- spec_df[-inTrain,]
  
  #tune model: 10-fold cross-validation repeated 3 times
  ctrl <- trainControl(
    method = "repeatedcv",
    number = 10,
    sampling = 'down',
    repeats = 3)
  
  #Fit model. Note max iterations set to 10000 to allow model convergence
  plsFit <- train(
    Location ~ .,
    data = training,
    maxit = 100000,
    method = "pls",
    trControl = ctrl,
    tuneLength = ncomp)
  
  #variable importance
  vip = varImp(plsFit)
  
  bg = assign(paste0('bg', i), vip$importance$`Bison Gulch`)
  bg.vip <- cbind(bg.vip, get('bg'))
  
  es = assign(paste0('es', i), vip$importance$`Eagle Summit`)
  es.vip <- cbind(es.vip, get('es'))
  
  md = assign(paste0('md', i), vip$importance$`Murphy Dome B`)
  md.vip <- cbind(md.vip, get('md'))
  
  tm = assign(paste0('tm', i), vip$importance$`Twelve Mile`)
  tm.vip <- cbind(tm.vip, get('tm'))
  
  wda = assign(paste0('wda', i), vip$importance$`Wickersham Dome A`)
  wda.vip <- cbind(wda.vip, get('wda'))
  
  wdb = assign(paste0('wdb', i), vip$importance$`Wickersham Dome B`)
  wdb.vip <- cbind(wdb.vip, get('wdb'))
  
  #accuracy objects for determining n components
  a = assign(paste0('a', i), as.matrix(plsFit$results$Accuracy))
  a.fit <- cbind(a.fit, get('a'))
  
  #test model using the testing data partition (30% of data)
  plsClasses <- predict(plsFit, newdata = testing)
  
  #confusion/classification matrix objects to assess accuracy 
  cm = confusionMatrix(data = plsClasses, as.factor(testing$Location))
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
abline(v = 18, col = 'blue')
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
rownames(cm.sd) <- c('BG', 'ES', 'MD', 'TM', 'WDA', 'WDB')
colnames(cm.sd) <- c('BG', 'ES', 'MD', 'TM', 'WDA', 'WDB')
write.csv(cm.sd, file = 'Figures/cm_final/6_scans/standard deviations/Location_sd.csv')

#format matrix for plotting
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('BG', 'ES', 'MD', 'TM', 'WDA', 'WDB')
colnames(cm.total) <- c('BG', 'ES', 'MD', 'TM', 'WDA', 'WDB')

#save confusion matrix
write.csv(cm.total, "Figures/cm_final/6_scans/Location.csv")

#plot confusion matrix
cols = colorRampPalette(c('#f5f5f5', '#fe9929'))

par(mar = c(1,2,2,1), oma = c(1,1,3,1))
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
mtext("Reference", side = 2, line = -2, cex = 2.5)
mtext("Prediction", side = 3, cex = 2.5, at = 3.5, line = 2)

################################################################################
#Variable importance
################################################################################
vip_to_spec = function(x){
  t.vip = t(x)
  colnames(t.vip) <- seq(400,2400, by = 1)
  s.vip = as_spectra(t.vip)
}
bg.vip = bg.vip[,-1]
bg.vip.spec = vip_to_spec(bg.vip)

es.vip = es.vip[,-1]
es.vip.spec = vip_to_spec(es.vip)

md.vip = md.vip[,-1]
md.vip.spec = vip_to_spec(md.vip)

tm.vip = tm.vip[,-1]
tm.vip.spec = vip_to_spec(tm.vip)

wda.vip = wda.vip[,-1]
wda.vip.spec = vip_to_spec(wda.vip)

wdb.vip = wdb.vip[,-1]
wdb.vip.spec = vip_to_spec(wdb.vip)

#plot
par(mfrow = c(3,2))
plot(mean(bg.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = NA, main = 'Bison Gulch')
plot_quantile(bg.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)

plot(mean(es.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = NA, main = 'Eagle Summit')
plot_quantile(es.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)

plot(mean(md.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = NA, main = 'Murphy Dome')
plot_quantile(md.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)


plot(mean(tm.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = NA, main = 'Twelve Mile')
plot_quantile(tm.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)

plot(mean(wda.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = 'Variable Importance',
     xlab = 'Wavelength (nm)', main = 'Wickersham Dome A')
plot_quantile(wda.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)

plot(mean(wdb.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = 'Wavelength (nm)', main = 'Wickersham Dome B')
plot_quantile(wdb.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)

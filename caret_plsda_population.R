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

#Code for new populations
#s.m = as_spectra(as.matrix(spec_all))
#meta(s.m) = read.csv('metadata_3.csv', stringsAsFactors = F)
#spec_all = s.m

#remove any NaN values - mostly pertains to populations
spec_all = spec_all[!meta(spec_all)$GenePop_ID == "NaN",]

spec_mat = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)


#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$GenePop_ID)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$GenePop_ID"] <- "GenePop_ID"


################################################################################
#Run PLSDA
################################################################################

#Set number of components to be used
ncomp = 45

#create vectors, lists, and matrices to store metrics and loadings
accuracy <- c()
kappa <- c()
a.fit <- matrix(nrow = ncomp)
cm.list <- list()
da_et.vip = matrix(nrow=2001)
da_wdb.vip = matrix(nrow=2001)
do_bg.vip = matrix(nrow=2001)
do_et.vip = matrix(nrow=2001)
do_md.vip = matrix(nrow=2001)
do_wd.vip = matrix(nrow=2001)

#start of PLSDA code
for(i in 1:100){
  
  #create data partition: 80% of data for training, 20% for testing
  inTrain <- caret::createDataPartition(
    y = spec_df$GenePop_ID,
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
    GenePop_ID ~ .,
    data = training,
    maxit = 100000,
    method = "pls",
    trControl = ctrl,
    tuneLength = ncomp)
  
  #variable importance
  vip = varImp(plsFit)
  
  da_et = assign(paste0('da_et', i), vip$importance$DA_et)
  da_et.vip <- cbind(da_et.vip, get('da_et'))
  
  da_wdb = assign(paste0('da_wdb', i), vip$importance$DA_wdb)
  da_wdb.vip <- cbind(da_wdb.vip, get('da_wdb'))
  
  do_bg = assign(paste0('do_bg', i), vip$importance$DO_bgc)
  do_bg.vip <- cbind(do_bg.vip, get('do_bg'))
  
  do_et = assign(paste0('do_et', i), vip$importance$DO_et)
  do_et.vip <- cbind(do_et.vip, get('do_et'))
  
  do_md = assign(paste0('do_md', i), vip$importance$DO_mdb)
  do_md.vip <- cbind(do_md.vip, get('do_md'))
  
  do_wd = assign(paste0('do_wd', i), vip$importance$DO_wd)
  do_wd.vip <- cbind(do_wd.vip, get('do_wd'))
  
  #accuracy objects for determining n components
  a = assign(paste0('a', i), as.matrix(plsFit$results$Accuracy))
  a.fit <- cbind(a.fit, get('a'))
  
  #test model using the testing data partition (30% of data)
  plsClasses <- predict(plsFit, newdata = testing)
  
  #confusion/classification matrix objects to assess accuracy 
  cm = confusionMatrix(data = plsClasses, as.factor(testing$GenePop_ID))
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
     xlab = 'Component', xlim = c(1,50), main = 'Accuracy for Species_ID', 
     ylim = c(0,1))
arrows(x, a.lower, x, a.higher,length=0.05, angle=90, code=3)
abline(v = 18, col = 'blue')
abline(h = max(a.lower), col = "Red")
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
rownames(cm.sd) <- c('ESTM', 'WDB', 'BG', 'ESTM', 'MD', 'WD')
colnames(cm.sd) <- c('ESTM', 'WDB', 'BG', 'ESTM', 'MD', 'WD')
write.csv(cm.sd, file = 'Figures/cm_final/6_scans/standard deviations/GenePop_ID_sd.csv')

#format matrix for plotting
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('ESTM', 'WDB', 'BG', 'ESTM', 'MD', 'WD')
colnames(cm.total) <- c('ESTM', 'WDB', 'BG', 'ESTM', 'MD', 'WD')

#save confusion matrix
write.csv(cm.total, "Figures/cm_final/6_scans/GenePop_ID.csv")


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
da_et.vip = da_et.vip[,-1]
da_et.vip.spec = vip_to_spec(da_et.vip)

da_wdb.vip = da_wdb.vip[,-1]
da_wdb.vip.spec = vip_to_spec(da_wdb.vip)

do_bg.vip = do_bg.vip[,-1]
do_bg.vip.spec = vip_to_spec(do_bg.vip)

do_et.vip = do_et.vip[,-1]
do_et.vip.spec = vip_to_spec(do_et.vip)

do_md.vip = do_md.vip[,-1]
do_md.vip.spec = vip_to_spec(do_md.vip)

do_wd.vip = do_wd.vip[,-1]
do_wd.vip.spec = vip_to_spec(do_wd.vip)

#plot
par(mfrow = c(2,3))
plot(mean(da_et.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = NA, main = 'DAK-ET')
plot_quantile(da_et.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)



plot(mean(da_wdb.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = NA, main = 'DAK-WDB')
plot_quantile(da_wdb.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)


plot(mean(do_bg.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = NA, main = 'DAJ-BG')
plot_quantile(do_bg.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)



plot(mean(do_et.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = 'Wavelength (nm)', main = 'DAJ-ET')
plot_quantile(do_et.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)


plot(mean(do_md.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA,
     xlab = 'Wavelength (nm)', main = 'DAJ-MD')
plot_quantile(do_md.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)


plot(mean(do_wd.vip.spec), lwd = 2, lty = 1, col = rgb(0,0,0,1), 
     cex.lab = 1.25, ylim = c(0, 100), ylab = NA, 
     xlab = 'Wavelength (nm)', main = 'DAJ-WD')
plot_quantile(do_wd.vip.spec, total_prob = 0.95, col = rgb(0,0,0, 0.25), 
              border = FALSE, add = TRUE)

################################################################################
#Loadings plot
################################################################################
#Convert data into 'spectral' data for fancy graphing
comp_to_spec = function(x){
  t.comp = t(x)
  colnames(t.comp) <- seq(400,2400, by = 10)
  s.comp = as_spectra(t.comp)
}

component1 = comp_to_spec(loadings.c1)
component2 = comp_to_spec(loadings.c2)
component3 = comp_to_spec(loadings.c3)

#plot loadings for first 3 components
dev.new(width = 6, height = 8, unit = 'in')
par(mar = c(5,4,1,1), oma = c(1,1,1,1), mfrow = c(1,1))
plot(mean(component1), lwd = 2, lty = 1, col = rgb(1,0,0,1), 
     cex.lab = 1.5, ylim = c(-.2, .15), ylab = "Loading Values", 
     xlab = "Wavelength (nm)")
plot_quantile(component1, total_prob = 0.95, col = rgb(1, 0, 0, 0.25), 
              border = FALSE, add = TRUE)
plot_regions(component1, regions = default_spec_regions(), add = TRUE)
plot(mean(component2), lwd = 1.5, lty = 1, col = rgb(0,0,1,1), add = TRUE)
plot_quantile(component2, total_prob = 0.95, col = rgb(0, 0, 1, 0.25), 
              border = FALSE, add = TRUE)
plot(mean(component3), lwd = 1.5, lty = 1, col = "darkgreen", add = TRUE)
plot_quantile(component3, total_prob = 0.95, col = rgb(0, .5, 0, 0.25),
              border = FALSE, add = TRUE)
abline(h = 0, lty = 2, lwd = 1.5)
legend('bottomright',inset = .02,
       legend=c("Component 1", "Component 2", 'Component 3'),
       col=c(rgb(1,0,0,1), rgb(0,0,1,1), rgb(0,1,0,1)), 
       lty=1, cex=0.8, bg ='white')

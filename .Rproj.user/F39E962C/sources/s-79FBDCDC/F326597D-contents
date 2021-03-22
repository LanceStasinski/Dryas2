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
#Data setup 
################################################################################

#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")


#remove any NaN values - mostly pertains to populations
spec_all = spec_all[!meta(spec_all)$Species_ID == "NaN",]

spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s

#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$Species_ID)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$Species_ID"] <- "Species_ID"


################################################################################
#Run PLSDA
################################################################################

#Set number of components to be used
ncomp = 22

#create vectors, lists, and matrices to store metrics and loadings
accuracy <- c()
kappa <- c()
k.fit <- matrix(nrow = ncomp)
cm.list <- list()
da.vip = matrix(nrow=201)
do.vip = matrix(nrow=201)
dx.vip = matrix(nrow=201)

#start of PLSDA code
for(i in 1:100){

#create data partition: 70% of data for training, 30% for testing
inTrain <- caret::createDataPartition(
  y = spec_df$Species_ID,
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
  Species_ID ~ .,
  data = training,
  maxit = 100000,
  method = "pls",
  trControl = ctrl,
  tuneLength = ncomp)

#variable importance
vip = varImp(plsFit)

da = assign(paste0('da', i), vip$importance$DA)
da.vip <- cbind(da.vip, get('da'))

do = assign(paste0('do', i), vip$importance$DO)
do.vip <- cbind(do.vip, get('do'))

dx = assign(paste0('dx', i), vip$importance$DX)
dx.vip <- cbind(dx.vip, get('dx'))

#kappa objects for determining n components
k = assign(paste0('k', i), as.matrix(plsFit$results$Kappa))
k.fit <- cbind(k.fit, get('k'))

#test model using the testing data partition (30% of data)
plsClasses <- predict(plsFit, newdata = testing)

#confusion/classification matrix objects to assess accuracy 
cm = confusionMatrix(data = plsClasses, as.factor(testing$Species_ID))
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
#Kappa values for choosing the optimal number of components to use
################################################################################

k.total = k.fit[,-1]
kavg = as.matrix(rowMeans(k.total))
ksd = as.matrix(rowSds(k.total))

klower = kavg - ksd
khigher = kavg + ksd

#Graph to visually choose optimal number of components
x = 1:60
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, kavg, type = 'p', pch = 16, cex = .75, ylab = 'Kappa', 
     xlab = 'Component', xlim = c(1,60), main = 'Kappa for Species_ID')
arrows(x, klower, x, khigher,length=0.05, angle=90, code=3)
abline(v = 29, col = 'blue')
abline(h = max(klower), col = "Red")
legend('bottomright', legend = c('Mean', 'Maximum kappa','Best component'), 
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
rownames(cm.sd) <- c('DA', 'DO', 'DX')
colnames(cm.sd) <- c('DA', 'DO', 'DX')
write.csv(cm.sd, file = 'Figures/cm_final/standard deviations/Species_ID_sd_upsample2_small2.csv')

#format matrix for plotting
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('DA', 'DO', 'DX')
colnames(cm.total) <- c('DA', 'DO', 'DX')

#save confusion matrix
write.csv(cm.total, "Figures/cm_final/cm_Species_ID_upsample2_small2.csv")

#species + Species_ID special code
cm.total = read.csv("Figures/cm_final/cm_Species_ID_upsample2_small2.csv", stringsAsFactors = T)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- cm.total[,1]
cm.total = cm.total[,-1]
cm.total = mapply(cm.total, FUN = as.numeric)
cm.total = matrix(data = cm.total, ncol = 12, nrow = 12)
rownames(cm.total) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'TM', 
                        'MD', 'WDA', 'WDB', 'ES', 'TM', 'WDB')
colnames(cm.total) <- c('ES', 'TM', 'WDB', 'BG', 'ES', 'TM', 
                        'MD', 'WDA', 'WDB', 'ES', 'TM', 'WDB')

#plot confusion matrix
cm.total = read.csv("Figures/cm_final/cm_Species_ID_upsample2_small2.csv", stringsAsFactors = F)
rownames(cm.total) = cm.total[,1]
cm.total = as.matrix(cm.total[,-1])


cols = colorRampPalette(c('#f5f5f5', '#fe9929'))

par(mar = c(1,2,2,1), oma = c(1,1,3,1))
corrplot::corrplot(cm.total,
         is.corr = F, 
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
  colnames(t.vip) <- seq(400,2400, by = 10)
  s.vip = as_spectra(t.vip)
}
do.vip = do.vip[,-1]
do.vip.spec = vip_to_spec(do.vip)

da.vip = da.vip[,-1]
da.vip.spec = vip_to_spec(da.vip)

dx.vip = dx.vip[,-1]
dx.vip.spec = vip_to_spec(dx.vip)

#plot
par(mfrow = c(3,1))

plot(mean(da.vip.spec), lwd = 1.5, lty = 1, col = '#00B0F6', ylim = c(0, 100),
     ylab = "Variable Importance", xlab = NA, cex.lab = 1.5)
plot_quantile(da.vip.spec, total_prob = 0.95, col = rgb(0, 0.69, 0.965, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)

plot(mean(do.vip.spec), lwd = 2, lty = 1, col = '#F8766D', 
     cex.lab = 1.5, ylim = c(0, 100), ylab = "Variable Importance", 
     xlab = NA)
plot_quantile(do.vip.spec, total_prob = 0.95, col = rgb(0.972549,0.4627451,0.427451, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)


plot(mean(dx.vip.spec), lwd = 1.5, lty = 1, col = rgb(0,0,0,1), ylim = c(0, 100),
     ylab = "Variable Importance", xlab = 'Wavelength (nm)', cex.lab = 1.5)
plot_quantile(dx.vip.spec, total_prob = 0.95, col = rgb(0, 0, 0, 0.25), 
              border = FALSE, add = TRUE)
abline(v = 1450, lty = 2, lwd = 1.5)
abline(v = 1940, lty = 2, lwd = 1.5)




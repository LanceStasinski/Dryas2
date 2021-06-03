#PLS with beta distribution for predicting ancestry
################################################################################
#Set up
################################################################################
library(caret)
library(plsRbeta)
library(pls)
library(spectrolab)
library(matrixStats)
library(tidyverse)
library(rlist)
library(parallel)
library(foreach)
library(doParallel)

################################################################################
#Set up
################################################################################
#spectra
spec_all= readRDS("Data/clean_all_6scans.rds")

#remove NAs
spec_all = spec_all[!meta(spec_all)$DA == "NaN",]

#prepare data for PLS
spectra.df = as.data.frame(spec_all)
spectra.m = as.matrix(spec_all)
spec_df = as.data.frame(spectra.m)
spec_df = cbind(spec_df, spectra.df$DA)
colnames(spec_df)[colnames(spec_df) == 'spectra.df$DA'] <- 'DA'

################################################################################
#Training and testing sets - all
################################################################################
numCores = detectCores()
registerDoParallel(numCores)

#parallelized plsBeta regression
plsFit = foreach (i = 1:5) %dopar% {
  library(plsRbeta)
  m = plsRbeta::PLS_beta_kfoldcv_formula(DA~., data = spec_df, nt = 60,
                           modele = 'pls-beta', K = 10, NK = 1,
                           verbose = T, random = T)
}
saveRDS(plsFit, 'Models/plsBeta/plsFit.rds')

#generate pls info statistics
pls_info_list = list()
for(i in 1:5) {
 pls_info = plsRbeta::kfolds2CVinfos_beta(plsFit[[i]])
 info = assign(paste0("pls_info", i), pls_info)
 pls_info_list = list.append(pls_info_list, get('info'))
}
saveRDS(plsFit.info, 'Models/plsBeta/pls_info_list.rds')

  


#Convert info to data frame, calculate RMSE from RSS
for(i in 1:5){
  info.df = as.data.frame(plsFit.info[[i]])
  info.df$RMSE <- sqrt(info.df$RSS_Y/460)
  assign(paste0('info.df_', i), info.df)
}

################################################################################
#Determine optimal number of components 
################################################################################
####################
#RMSE - best component probably 77
####################

#create RMSE matrix
rmse_total = Reduce(cbind, list(info.df_1$RMSE, info.df_2$RMSE, info.df_3$RMSE,
                                info.df_4$RMSE, info.df_5$RMSE))
rownames(rmse_total) <- c(0:80) #number of components
colnames(rmse_total) <- c(1:5) #number of repeats

#Calculate mean and standard deviation by component
mean_rmse = as.matrix(rowMeans(rmse_total))
sd_rmse = as.matrix(rowSds(rmse_total))

#Upper and lower bounds for plot
rmse_lower = mean_rmse - sd_rmse
rmse_upper = mean_rmse + sd_rmse

#Graph to visually choose optimal number of components
x = 0:80
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, mean_rmse, type = 'p', pch = 16, cex = .75, ylab = 'RMSE', 
     xlab = 'Component', xlim = c(0,80), main = 'RMSE vs Component')
arrows(x, rmse_lower, x, rmse_upper,length=0.05, angle=90, code=3)
abline(v = 77, col = 'blue')
abline(h = min(rmse_upper), col = "Red")
legend('topright', legend = c('Mean', 'Lowest RMSE','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))

####################
#AIC - best probably component 54 (first component within 2 AIC of lowest AIC)
####################
#create aic matrix
aic_total = Reduce(cbind, list(info.df_1$AIC, info.df_2$AIC, info.df_3$AIC,
                                info.df_4$AIC, info.df_5$AIC))
rownames(aic_total) <- c(0:80) #number of components
colnames(aic_total) <- c(1:5) #number of repeats

#Calculate mean and standard deviation by component
mean_aic = as.matrix(rowMeans(aic_total))
sd_aic = as.matrix(rowSds(aic_total))

#lowest AIC
best_aic = min(mean_aic)

#Upper and lower bounds for plot
aic_lower = mean_aic - sd_aic
aic_upper = mean_aic + sd_aic

#Graph to visually choose optimal number of components
x = 0:80
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, mean_aic, type = 'p', pch = 16, cex = .75, ylab = 'AIC', 
     xlab = 'Component', xlim = c(0,80), main = 'AIC vs Component')
arrows(x, aic_lower, x, aic_upper,length=0.05, angle=90, code=3)
abline(v = 54, col = 'blue')
abline(h = min(mean_aic) + 2, col = "Red")
legend('topright', legend = c('Mean', 'Lowest AIC + 2','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))

####################
#BIC - best probably component 46 (first component within 2 BIC of lowest BIC)
####################
#create BIC matrix
BIC_total = Reduce(cbind, list(info.df_1$BIC, info.df_2$BIC, info.df_3$BIC,
                               info.df_4$BIC, info.df_5$BIC))
rownames(BIC_total) <- c(0:80) #number of components
colnames(BIC_total) <- c(1:5) #number of repeats

#Calculate mean and standard deviation by component
mean_BIC = as.matrix(rowMeans(BIC_total))
sd_BIC = as.matrix(rowSds(BIC_total))

#lowest BIC
best_BIC = min(mean_BIC)

#Upper and lower bounds for plot
BIC_lower = mean_BIC - sd_BIC
BIC_upper = mean_BIC + sd_BIC

#Graph to visually choose optimal number of components
x = 0:80
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, mean_BIC, type = 'p', pch = 16, cex = .75, ylab = 'BIC', 
     xlab = 'Component', xlim = c(0,80), main = 'BIC vs Component')
arrows(x, BIC_lower, x, BIC_upper,length=0.05, angle=90, code=3)
abline(v = 46, col = 'blue')
abline(h = min(mean_BIC) + 2, col = "Red")
legend('topright', legend = c('Mean', 'Lowest BIC + 2','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))

################################################################################
#Create plot of predicted ancestry vs actual ancestry
################################################################################
#obtain predicted ancestry for each repeat and group from model object

for(NK in 1:5){
  results.list = list()
  for(k in 1:10){
  matr = as.matrix(plsFit$results_kfolds[[NK]][[k]][,54])
  result.matrix = assign(paste0('matr',k), matr)
  results.list <- list.append(results.list, get('result.matrix'))
  }
  full.mat = Reduce(rbind, results.list)
  full.mat = full.mat[row.names(spec_df),,drop = F]
  assign(paste0('full.mat_', NK), full.mat)
}

#average the predictions
mean.predictions = as.data.frame((full.mat_1 + full.mat_2 + full.mat_3 +
                                  full.mat_4 + full.mat_5)/5)

#Add pertinent metadata
colnames(mean.predictions) <- c("Predictions")
mean.predictions$DA <- spectra.df$DA
mean.predictions$Species_ID <- spectra.df$Species_ID
mean.predictions$Name <- spectra.df$Name


#find median and range of predicted values
l = list(mean.predictions$Name)
med = aggregate(mean.predictions$Predictions, by = l, median)
DA = aggregate(DA~Name, data = mean.predictions, mean)
DA$Predictions = med[,2]
max = aggregate(mean.predictions$Predictions, by = l, max)
DA$max = max[,2]
min = aggregate(mean.predictions$Predictions, by = l, min)
DA$min = min[,2]

#Add colors
DA$color = 'black'
DA$color[DA$DA > .7] = "#00B0F6"
DA$color[DA$DA < .3] = "#F8766D"

#plot actual ancestry vs predicted ancestry
dev.new(width = 6, height = 6, unit = 'in')
par(mfrow = c(1,1))
plot(DA$DA, DA$Predictions, cex.lab = 1.5,
     xlab = "Actual Dryas alaskensis ancestry",
     ylab = "Predicted Dryas alaskensi ancestry", pch = 16)
lines(x = c(-2, 2), y = c(-2, 2), lty=2)
arrows(DA$DA, DA$min, DA$DA, DA$max,length=0.05, angle=90, code=0, col = DA$color)
points(DA$DA, DA$Predictions,
       col = DA$color, pch = 16)
legend("bottomright", inset = 0.01,
       legend=c("D. ajanensis", "Hybrid", "D. alaskensis", "1:1 line"), 
       text.font = c(3,1,3,1),
       col=c("#F8766D", "black", "#00B0F6", "black"), pch = c(16,16,16,NA),
       lty = c(NA,NA,NA,2), bg = 'white')
#extra
c.rmse = caret::RMSE(mean.predictions$Predictions, mean.predictions$DA)
c.r2 = caret::R2(mean.predictions$Predictions, mean.predictions$DA)

a = caret::RMSE(full.mat_1, spectra.df$DA)
b = caret::RMSE(full.mat_2, spectra.df$DA)
c = caret::RMSE(full.mat_3, spectra.df$DA)
d = caret::RMSE(full.mat_4, spectra.df$DA)
e = caret::RMSE(full.mat_5, spectra.df$DA)
rmse = c(a,b,c,d,e)
sd(rmse)

f = caret::R2(full.mat_1, spectra.df$DA)
g = caret::R2(full.mat_2, spectra.df$DA)
h = caret::R2(full.mat_3, spectra.df$DA)
i = caret::R2(full.mat_4, spectra.df$DA)
j = caret::R2(full.mat_5, spectra.df$DA)
r2 = c(f,g,h,i,j)
sd(r2)

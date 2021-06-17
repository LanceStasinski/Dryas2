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
spec_all = resample(spec_all, seq(400, 2400, 2))

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
  m = plsRbeta::PLS_beta_kfoldcv_formula(DA~., data = spec_df, nt = 50,
                           modele = 'pls-beta', K = 10, NK = 1,
                           verbose = T, random = T)
}
saveRDS(plsFit, 'Models/plsBeta/plsFit_2nm.rds')

plsFit = readRDS('Models/plsBeta/plsFit.rds')

pls_info = plsRbeta::kfolds2CVinfos_beta(plsFit[[1]])
saveRDS(pls_info, 'Models/plsBeta/pls_info.rds')
pls_info = readRDS('Models/plsBeta/pls_info.rds')

#generate pls info statistics
pls_info_list = list()
for(i in 1:5) {
 pls_info = plsRbeta::kfolds2CVinfos_beta(plsFit[[i]])
 info = assign(paste0("pls_info", i), pls_info)
 pls_info_list = list.append(pls_info_list, get('info'))
}
saveRDS(plsFit.info, 'Models/plsBeta/pls_info_list.rds')

  


#Convert info to data frame, calculate RMSE from RSS

info.df = as.data.frame(pls_info)
info.df$RMSE <- sqrt(info.df$RSS_Y/460)


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

for(NK in 1:5) {
  results.list = list()
  for(k in 1:10){
    if (ncol(plsFit[[NK]]$results_kfolds[[1]][[k]]) < 50) next
      
    matr = as.matrix(plsFit[[NK]]$results_kfolds[[1]][[k]][,50])
    result.matrix = assign(paste0('matr',k), matr)
    results.list <- list.append(results.list, get('result.matrix'))
  }
  full.mat = Reduce(rbind, results.list)
  assign(paste0('full.mat_', NK), full.mat)
}


#finding the average prediction really only works when errors do not occur in 
#PLS that result in termination of the iteration before it reaches the full 
#component cap
#average the predictions
mean.predictions = as.data.frame((full.mat_1 + full.mat_2 + full.mat_3 +
                                  full.mat_4 + full.mat_5)/5)

#plot mean and standard deviation of one repetition 
mean.predictions = full.mat_4

remove_spec = setdiff(spectra.df$sample_name, rownames(mean.predictions))
spec_df2 = spectra.df[!spectra.df$sample_name %in% remove_spec, ]

mean.predictions = as.data.frame(mean.predictions)
colnames(mean.predictions) <- c("Predictions")
mean.predictions = mean.predictions[spec_df2$sample_name,, drop = F]

mean.predictions$DA <- spec_df2$DA
mean.predictions$Species_ID <- spec_df2$Species_ID
mean.predictions$Name <- spec_df2$Name

l = list(mean.predictions$Name)
pred_mean = aggregate(mean.predictions$Predictions, by = l, FUN = mean)
rownames(pred_mean) = pred_mean[,1]
pred_mean = pred_mean[,-1]

pred_sd = aggregate(mean.predictions$Predictions, by = l, FUN = sd)
pred_sd2 = pred_sd[,-1]
pred_stats = cbind(pred_mean, pred_sd2)
colnames(pred_stats) = c('mean', 'sd')
rownames(pred_stats) = pred_sd[,1]
pred_stats = as.data.frame(pred_stats)
pred_stats$lower = pred_stats$mean - pred_stats$sd
pred_stats$higher = pred_stats$mean + pred_stats$sd

DA = aggregate(DA~Name, data = mean.predictions, mean)
DA = DA[,-1]
pred_stats$DA = DA

#plot
#colors
pred_stats$color = 'black'
pred_stats$color[pred_stats$DA > .7] = "#00B0F6"
pred_stats$color[pred_stats$DA < .3] = "#F8766D"

dev.new(width = 6, height = 6, unit = 'in')
par(mfrow = c(1,1))
plot(pred_stats$DA, pred_stats$mean,
     cex.lab = 1.5,
     xlab = "Actual Dryas alaskensis ancestry",
     ylab = "Predicted Dryas alaskensis ancestry",
     pch = 16,
     ylim = c(-0.1, 1.1),
     xlim = c(-0.1, 1.1))
lines(x = c(-2,2), y = c(-2,2), lty = 2)
arrows(pred_stats$DA, pred_stats$lower,
       pred_stats$DA, pred_stats$higher,
       angle = 90,
       length = 0.05,
       code = 3,
       col = pred_stats$color)
points(pred_stats$DA, pred_stats$mean,
       col = pred_stats$color,
       pch = 16)
legend("bottomright", inset = 0.01,
       legend=c("D. ajanensis", "Hybrid", "D. alaskensis", "1:1 line"), 
       text.font = c(3,1,3,1),
       col=c("#F8766D", "black", "#00B0F6", "black"), pch = c(16,16,16,NA),
       lty = c(NA,NA,NA,2), bg = 'white')


#plot median prediction and full range of predictions
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
plot(DA$DA, DA$Predictions,
     cex.lab = 1.5,
     xlab = "Actual Dryas alaskensis ancestry",
     ylab = "Predicted Dryas alaskensis ancestry",
     pch = 16)
lines(x = c(-2, 2), y = c(-2, 2), lty=2)
arrows(DA$DA, DA$min,
       DA$DA, DA$max,
       length=0.05,
       angle=90, 
       code=0, 
       col = DA$color)
points(DA$DA, DA$Predictions,
       col = DA$color, 
       pch = 16)
legend("bottomright", inset = 0.01,
       legend=c("D. ajanensis", "Hybrid", "D. alaskensis", "1:1 line"), 
       text.font = c(3,1,3,1),
       col=c("#F8766D", "black", "#00B0F6", "black"), pch = c(16,16,16,NA),
       lty = c(NA,NA,NA,2), bg = 'white', cex = .75)
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

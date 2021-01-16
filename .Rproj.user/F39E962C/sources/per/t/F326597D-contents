#Caret PLSDA
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
#Data setup !!!!NOTE: use ctrl+f to find a replace the field to be classified!!!
################################################################################

#data
spec_all = readRDS("Clean-up/Clean_spectra/clean_all.rds")
spec_all = spec_all[!meta(spec_all)$Species_ID == 'DX',]

#Code for new populations
s.m = as_spectra(as.matrix(spec_all))
meta(s.m) = read.csv('metadata_2.csv', stringsAsFactors = F)
spec_all = s.m

#remove any NaN values - mostly pertains to populations
spec_all = spec_all[!meta(spec_all)$sp_loc == "NaN",]

spec_all.m = as.matrix(spec_all)
spec_all.df = as.data.frame(spec_all)

#Resample by every 10 nm
spec_small = resample(spec_all, seq(400, 2400, by = 10))
spec_mat_s = as.matrix(spec_small)
spec_mat = spec_mat_s

#combine relavant meta data to matrix
spec_df = as.data.frame(spec_mat)
spec_df = cbind(spec_df, spec_all.df$sp_loc)
colnames(spec_df)[colnames(spec_df) == "spec_all.df$sp_loc"] <- "sp_loc"

################################################################################
#Run PLSDA
################################################################################

#Set number of components to be used
ncomp = 33

#create vectors, lists, and matrices to store metrics and loadings
accuracy <- c()
kappa <- c()
k.fit <- matrix(nrow = ncomp)
cm.list <- list()
loadings.c1 <- matrix(nrow = 201)
loadings.c2 <- matrix(nrow = 201)
loadings.c3 <- matrix(nrow = 201)

#start of PLSDA code
for(i in 1:100){

#create data partition: 70% of data for training, 30% for testing
inTrain <- caret::createDataPartition(
  y = spec_df$sp_loc,
  p = .7,
  list = FALSE
)

training <- spec_df[inTrain,]
testing <- spec_df[-inTrain,]

#tune model: 10-fold cross-validation repeated 3 times
ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3)

#Fit model. Note max iterations set to 10000 to allow model convergence
plsFit <- train(
  sp_loc ~ .,
  data = training,
  maxit = 10000,
  method = "pls",
  trControl = ctrl,
  tuneLength = ncomp)

#kappa objects for determining n components
k = assign(paste0('k', i), as.matrix(plsFit$results$Kappa))
k.fit <- cbind(k.fit, get('k'))

#loadings
c1 = assign(paste0('c1',i), as.matrix(plsFit$finalModel$loadings[,1]))
loadings.c1 <- cbind(loadings.c1, get('c1'))

c2 = assign(paste0('c2',i), as.matrix(plsFit$finalModel$loadings[,2]))
loadings.c2 <- cbind(loadings.c2, get('c2'))

c3 = assign(paste0('c3',i), as.matrix(plsFit$finalModel$loadings[,3]))
loadings.c3 <- cbind(loadings.c3, get('c3'))

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
#Kappa values for choosing the optimal number of components to use
################################################################################

k.total = k.fit[,-1]
kavg = as.matrix(rowMeans(k.total))
ksd = as.matrix(rowSds(k.total))

klower = kavg - ksd
khigher = kavg + ksd

#Graph to visually choose optimal number of components
x = 1:50
par(mar = c(5.1, 4.1, 4.1, 2.1), oma = c(5.1, 4.1, 4.1, 2.1))
plot(x, kavg, type = 'p', pch = 16, cex = .75, ylab = 'Kappa', 
     xlab = 'Component', xlim = c(1,60), main = 'Kappa for sp_loc')
arrows(x, klower, x, khigher,length=0.05, angle=90, code=3)
abline(v = 38, col = 'blue')
abline(h = max(klower), col = "Red")
legend('bottomright', legend = c('Mean', 'Maximum kappa','Best component'), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), col = c('black', 'red', 'blue'))

################################################################################
#Confusion/Classification Matrices
################################################################################
#take average of 100 confusion matrices, reorient matrix
cm.total = Reduce('+', cm.list)/100
cm.total = t(cm.total)
cm.total = cm.total/rowSums(cm.total)

#format matrix for plotting
cm.total = as.data.frame(cm.total)
cm.total = cm.total %>% replace_with_na_all(condition = ~.x == 0)
cm.total = as.matrix(cm.total)
rownames(cm.total) <- c('ES-TM', 'WDB', 'BG', 'ES-TM', 'MD', 'WD')
colnames(cm.total) <- c('ES-TM', 'WDB', 'BG', 'ES-TM', 'MD', 'WD')

#save confusion matrix
write.csv(cm.total, "Figures/cm_final/cm_sp_loc_mean_newpops.csv")

#species + sp_loc special code
cm.total = read.csv("Figures/cm_final/cm_sp_loc_mean.csv", stringsAsFactors = T)
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
cols = colorRampPalette(c('#f5f5f5', '#b35806'))

par(mar = c(1,2,4,1), oma = c(1,1,3,1))
corrplot(cm.total,
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
mtext("Reference", side = 2, line = 2, cex = 2.5)
mtext("Prediction", side = 3, cex = 2.5, at = 3.5, line = 6)


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
